import torch
import torch.nn as nn
import dgl
from dgl.nn.pytorch.glob import AvgPooling, SumPooling
from torch.cuda.amp import autocast
import logging
import os
from mlmm.utils.utils import check_value
from mlmm.nn.blocks import MLP
if "use_fp16" in os.environ:
    if os.environ['use_fp16'] == 'True':
        use_fp16=True
    elif os.environ['use_fp16'] == 'False':
        use_fp16=False
    else:
        raise AssertionError("wrong setting, use_fp16 can only be True or False")
else:
    use_fp16=False

__all__=["graph_norm","gen_read_out"]

class graph_norm(nn.Module):
    def __init__(self,in_feature):
        super(graph_norm,self).__init__()
        self.alpha=torch.nn.Parameter(torch.ones(in_feature))
        self.beta=torch.nn.Parameter(torch.zeros(in_feature))
        self.gamma=torch.nn.Parameter(torch.ones(in_feature))
        self.pooling = AvgPooling()
        if use_fp16:
            self.epsilon = 1e-6
        else:
            self.epsilon = 1e-6
        self.g = None

    def forward(self,g,feat):
        logger = logging.getLogger('rl.'+__name__)
        self.g = g.local_var()

        if len(feat.shape)==3:
            feat = feat.reshape(feat.shape[0],-1)
            dipole=True
        else:
            dipole=False

        u = self.pooling(g,feat)
        logger.debug(f"u:{u.sum()}")
        u2 = dgl.broadcast_nodes(self.g,torch.matmul(u,torch.diag(self.alpha)))
        logger.debug(f"u2:{u2.sum()}")
        sigma = self.pooling(self.g,(feat - u2)**2) 
        logger.debug(f"sigma1:{sigma.sum()}")
        if check_value(sigma):
            raise Exception('error')
        sigma = dgl.broadcast_nodes(self.g,torch.pow(sigma + self.epsilon,0.5))
        logger.debug(f"sigma2:{sigma.sum()}")
        y = feat - u2
        y = y/(sigma + self.epsilon)
        y=y*self.gamma[None]
        y=y+self.beta[None]

        logger.debug(f"graph norm: {y.sum()}")
        if check_value(y):
            logger.critical('value error: nan or inf')
            assert check_value(y)==0 
        if dipole:
            y = y.reshape(feat.shape[0],-1,3)
        return y

class graph_norm_v(nn.Module):
    def __init__(self,in_feature, activation=None):
        super(graph_norm_v,self).__init__()
        self.alpha=torch.nn.Parameter(torch.ones(in_feature))
        self.beta=torch.nn.Parameter(torch.zeros(in_feature))
        self.gamma=torch.nn.Parameter(torch.ones(in_feature))
        self.zeta =torch.nn.Parameter(torch.ones(in_feature))

        self.mlp = MLP(in_feature,in_feature,activation=activation)

        self.pooling = AvgPooling()
        if use_fp16:
            self.epsilon = 1e-6
        else:
            self.epsilon = 1e-6
        self.g = None

    def forward(self,g,feat_v):
        logger = logging.getLogger('rl.'+__name__)
        self.g = g.local_var()

        feat = torch.norm(feat_v+1e-8,dim=-1)

        u = self.pooling(g,feat)
        logger.debug(f"u:{u.sum()}")
        u2 = dgl.broadcast_nodes(self.g,torch.matmul(u,torch.diag(self.alpha)))
        logger.debug(f"u2:{u2.sum()}")
        sigma = self.pooling(self.g,(feat - u2)**2) 
        logger.debug(f"sigma1:{sigma.sum()}")
        if check_value(sigma):
            raise Exception('error')
        sigma = dgl.broadcast_nodes(self.g,torch.pow(sigma + self.epsilon,0.5))
        logger.debug(f"sigma2:{sigma.sum()}")
        y = feat - u2
        y = y/(sigma + self.epsilon)
        y=y*self.gamma[None]
        y=y+self.beta[None]

        feat_v = self.mlp(y)[...,None]*(feat_v/(feat+self.zeta[None]+1e-8)[...,None])

        logger.debug(f"graph norm: {y.sum()}")
        if check_value(y):
            logger.critical('value error: nan or inf')
            assert check_value(y)==0 
        return feat_v

    
class gen_read_out(nn.Module):
    def __init__(self,agg='softmax'):
        super(gen_read_out,self).__init__()
        
        if agg=='powermean':
            self.p=torch.nn.Parameter(torch.Tensor([1]))
        elif agg=='softmax':
            self.p=torch.nn.Parameter(torch.Tensor([0]))
        self.beta=torch.nn.Parameter(torch.Tensor([1]))
        self.pooling = SumPooling()
        self.agg = agg
        if use_fp16:
            self.eps = 1e-7
        else:
            self.eps = 1e-25
    
    @autocast(enabled=use_fp16) 
    def forward(self,g,feature):
        with g.local_scope():
            num_nodes = g.number_of_nodes()/g.batch_size
            if self.agg=='softmax':
                pre=num_nodes/(1+self.beta*(num_nodes-1)+self.eps)
                g.ndata['h2'] = feature*self.p
                coef = dgl.softmax_nodes(g,'h2')
                out = self.pooling(g,coef*feature)
                return pre*out 

            elif self.agg=='powermean':
                pre=1/(1+self.beta*(num_nodes-1)+self.eps)
                out=torch.pow(pre*self.pooling(g,torch.pow(feature,self.p)),1/self.p)
                return out 



    
