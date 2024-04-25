import torch
from mlmm.nn.blocks import MLP, Dense
from dgl.nn.pytorch.glob import SumPooling
from typing import Union
from mlmm.nn.activations import *
from mlmm.utils.utils import add_sn

__all__ = ['ReadOut_ML', 'ReadOut_ML2']


class ReadOut_ML(torch.nn.Module): 
    """
    ReadOut function based on MLP
     Args:
       n_atom_basis(int): number of features of atom
       outdim_local(int): number of local properties predicted
       outdim_global(int): number of global properties predicted
       activation(str): activation function
       n_layers(int): the number of Fully connected layers
       n_hidden(list): hidden features in intermidate layers
       pooling(nn.Module): pooling method
    """
    def __init__(self,
                 n_atom_basis: int = 32,
                 m_max_order: int = 0,
                 outdim_global: int = 1,
                 activation: str = 'swish',
                 n_layers: int = 2,
                 n_hidden: Union[list,None] = None,
                 pooling: str = 'sumpooling',
                 dropout: bool = False,
                 spectral_norm: bool = False,
                 p: float = 0.1):
        super().__init__()
        
        self.name = 'mlp'

        if pooling == 'sumpooling':
            self.pooling = SumPooling()

        self.out = torch.nn.ModuleDict()
        if activation == "swish":
            activation = swish
        elif activation == "shifted_softplus":
            activation = shifted_softplus
        
        if outdim_global != None:
            self.out['global'] = MLP(n_atom_basis,
                                     outdim_global,
                                     n_hidden=n_hidden,
                                     n_layers=n_layers,
                                     activation = activation,
                                     dropout=dropout,
                                     p=p)
            
            if spectral_norm:
                self.out['global'] = add_sn(self.out['global'])

        self.multipole_net = multipole(n_atom_basis=n_atom_basis,
                                       n_hidden=n_hidden,
                                       n_layers=n_layers,
                                       activation=activation,
                                       dropout=dropout,
                                       max_order=m_max_order,
                                       p = p
                                       )

        # if outdim_local != None:
            # self.out['local'] = MLP(n_atom_basis,
                                    # outdim_local, 
                                    # n_hidden=n_hidden, 
                                    # n_layers=n_layers,
                                    # activation=activation,
                                    # dropout=dropout,
                                    # p=p)

            # if spectral_norm:
                # self.out['local'] = add_sn(self.out['local'])


    def forward(self, g, feat, feat_v, graph_rep=None):
        # if 'local' in self.out:
            # local_prop = self.out['local'](feat)
        # else:
            # local_prop = None
        
        multipole_prop = self.multipole_net(feat, feat_v)

        if 'global' in self.out:
            if graph_rep is None:
               global_prop =  self.out['global'](feat)
               global_prop = self.pooling(g, global_prop)
            else:
               global_prop =  self.out['global'](graph_rep)
        else:
            global_prop = None
        
        return multipole_prop, global_prop


class multipole(torch.nn.Module):
    def __init__(self,n_atom_basis,n_hidden,n_layers,activation, dropout, p, max_order = 1):
        super().__init__()
        self.out = torch.nn.ModuleDict()
        self.max_order = max_order

        if max_order == 0:
            outdim = 1
        else:
            outdim = torch.arange(self.max_order+1).sum() + self.max_order,

        self.out['s'] = MLP(n_atom_basis,
                            outdim,
                            n_hidden=n_hidden, 
                            n_layers=n_layers,
                            activation=activation,
                            dropout=dropout,
                            p=p)
                               

        for i in range(max_order):
            self.out['M'+str(i+1)] = Dense(n_atom_basis,i+1,bias=False)

    def forward(self, feat, feat_v):
        out = {}
        s = self.out['s'](feat)
        out['M0'] = s[:,0]

        idx = 1
        decomp = {}
        for i in range(self.max_order):
            decomp['M'+str(i+1)] = s[:,idx:idx+i+1][...,None]*self.out['M'+str(i+1)](feat_v.permute(0,2,1)).permute(0,2,1)
            idx = idx + 1 + i + 1 

        if self.max_order > 0 :
            out['M1'] = decomp['M1']
        
        if self.max_order > 1:
            if self.max_order == 2:
                out2 = tensor_product_multipole(decomp,s[:,2:3])
            elif self.max_order == 3:
                out2 = tensor_product_multipole(decomp,s[:,[2,5]],max_order = 3)
            out.update(out2)
        return out

def tensor_product_multipole(decomp, s, max_order = 2):

    out = {}
    I = torch.eye(3,dtype=torch.float32,device=s.device)
    out['M2'] = I[None]*s[:,0][:,None,None] + \
                s[:,1][:,None,None]*(decomp['M2'][:,0,:,None]*decomp['M2'][:,1,None,:] + \
                                     decomp['M2'][:,0,None,:]*decomp['M2'][:,1,:,None])

    if max_order > 2:
        out['M3'] = I[None,None,:,:]*I[None,:,None,:]*I[None,:,:,None]*s[:,2][:,None,None,None] + \
                    s[:,3][:,None,None,None]*(I[None,None,:,:]*decomp['M3'][:,2,:,None,None] + \
                                              I[None,:,None,:]*decomp['M3'][:,2,None,:,None] + \
                                              I[None,:,:,None]*decomp['M3'][:,2,None,None,:]) + \
                    s[:,4][:,None,None,None]*(I[None,None,:,:]*decomp['M3'][:,0,:,None,None] + \
                                              I[None,:,None,:]*decomp['M3'][:,0,None,:,None] + \
                                              I[None,:,:,None]*decomp['M3'][:,0,None,None,:]) + \
                    s[:,5][:,None,None,None]*(I[None,None,:,:]*decomp['M3'][:,1,:,None,None] + \
                                              I[None,:,None,:]*decomp['M3'][:,1,None,:,None] + \
                                              I[None,:,:,None]*decomp['M3'][:,1,None,None,:]) + \
                    s[:,6][:,None,None,None]*(decomp['M3'][:,0,:,None,None]*decomp['M3'][:,1,None,:,None]*decomp['M3'][:,2,None,None,:]+ \
                                              decomp['M3'][:,0,:,None,None]*decomp['M3'][:,1,None,None,:]*decomp['M3'][:,2,None,:,None]+ \
                                              decomp['M3'][:,0,None,:,None]*decomp['M3'][:,1,:,None,None]*decomp['M3'][:,2,None,None,:]+ \
                                              decomp['M3'][:,0,None,:,None]*decomp['M3'][:,1,None,None,:]*decomp['M3'][:,2,:,None,None]+ \
                                              decomp['M3'][:,0,None,None,:]*decomp['M3'][:,1,:,None,None]*decomp['M3'][:,2,None,:,None]+ \
                                              decomp['M3'][:,0,None,None,:]*decomp['M3'][:,1,None,:,None]*decomp['M3'][:,2,:,None,None])
    return out

class ReadOut_ML2(torch.nn.Module): 
    """
    ReadOut function based on MLP
     Args:
       n_atom_basis(int): number of features of atom
       outdim_local(int): number of local properties predicted
       outdim_global(int): number of global properties predicted
       activation(str): activation function
       n_layers(int): the number of Fully connected layers
       n_hidden(list): hidden features in intermidate layers
       pooling(nn.Module): pooling method
    """
    def __init__(self,
                 n_atom_basis: int = 32,
                 m_max_order: int = 0,
                 outdim_s: int = 1,
                 outdim_v: int = 0,
                 activation: str = 'swish',
                 n_layers: int = 2,
                 n_hidden: Union[list,None] = None,
                 pooling: str = 'sumpooling',
                 dropout: bool = False,
                 spectral_norm: bool = False,
                 p: float = 0.1):
        super().__init__()
        
        self.name = 'mlp'

        if pooling == 'sumpooling':
            self.pooling = SumPooling()

        if activation == "swish":
            activation = swish
        elif activation == "shifted_softplus":
            activation = shifted_softplus

        self.outdim_s = outdim_s
        self.outdim_v = outdim_v
        self.m_max_order = m_max_order

        outdim_s = outdim_s + 1
        outdim_v_ = outdim_v

        outdim_v_ = outdim_v_ + int((torch.arange(m_max_order)+1).sum())
        
        self.out_s1 = MLP(n_atom_basis,
                          outdim_s,
                          n_hidden=n_hidden,
                          n_layers=n_layers,
                          activation = activation,
                          dropout=dropout,
                          p=p)

        if m_max_order>=1:
            if m_max_order == 1:
                outdim_s_ = 1
            elif m_max_order == 2:
                outdim_s_ = 3
            elif m_max_order == 3:
                outdim_s_ = 8
            self.out_s2 = MLP(n_atom_basis*2,
                              outdim_s_,
                              n_hidden=n_hidden,
                              n_layers=n_layers,
                              activation = activation,
                              dropout=dropout,
                              p=p)
        
        
        if spectral_norm:
            self.out_s1 = add_sn(self.out_s1)
            if m_max_order > 0:
                self.out_s2 = add_sn(self.out_s2)

        if outdim_v_ > 0:
            self.out_v = Dense(n_atom_basis,outdim_v_,bias=False)

    def forward(self, g, feat, feat_v, graph_rep=None):
        
        s = self.out_s1(feat)
        if self.m_max_order > 0:
            s2 = self.out_s2(torch.cat([feat,torch.norm(feat_v,dim=-1)],dim=-1))

        multipole_prop = {}
        multipole_prop['M0'] = s[:,self.outdim_s]

        if self.m_max_order>0:
            v = self.out_v(feat_v.permute(0,2,1)).permute(0,2,1)
            idx = 2
            idx_v = 0
            decomp = {}
            for i in range(self.m_max_order):
                decomp['M'+str(i+1)] = v[:,self.outdim_v + idx_v:self.outdim_v + idx_v + i + 1]
                idx = idx + 1 + i + 1 
                idx_v = idx_v + i + 1

        if self.m_max_order > 0:
            multipole_prop['M1'] = s2[:,0][...,None]*v[:,self.outdim_v]
            if self.m_max_order == 2:
                out = tensor_product_multipole(decomp,s2[:,1:3],max_order=self.m_max_order)
            elif self.m_max_order == 3:
                out = tensor_product_multipole(decomp,s2[:,1:],max_order=self.m_max_order)
            multipole_prop.update(out)

        global_prop = self.pooling(g, s[:,:self.outdim_s])

        return multipole_prop, global_prop

