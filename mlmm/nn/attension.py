from mlmm.nn.base import Dense
import torch
import torch.nn as nn
from mlmm.nn.activations import shifted_softplus, swish
import logging
from dgl.ops import edge_softmax
import dgl.function as fn
from torch.cuda.amp import autocast
from mlmm.utils.utils import check_value
import os
if "use_fp16" in os.environ:
    if os.environ['use_fp16'] == 'True':
        use_fp16=True
    elif os.environ['use_fp16'] == 'False':
        use_fp16=False
    else:
        raise AssertionError("wrong setting, use_fp16 can only be True or False")
else:
    use_fp16=False

class attension(nn.Module):
    def __init__(self,n_atom_basis,n_attension_basis):
        super(attension, self).__init__()
        self.dense1 = Dense(n_atom_basis, n_attension_basis, bias=True)
        self.dense2 = Dense(n_atom_basis, n_attension_basis, bias=True)
        self.dense3 = Dense(n_atom_basis,1,bias=True)

    @autocast(enabled=use_fp16)    
    def forward(self,g,feat):

        with g.local_scope():
            g.ndata['att1']=self.dense1(feat)
            g.ndata['att2']=self.dense2(feat)
            g.apply_edges(fn.u_add_v('att1','att2','alpha'))
            return edge_softmax(g,swish(g.edata['alpha']),norm_by='dst')
        

        
