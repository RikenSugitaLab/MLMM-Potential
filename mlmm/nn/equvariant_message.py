import torch
from mlmm.nn.blocks import MLP
from mlmm.nn.activations import swish
from mlmm.nn.update_coor import Update_coor
from mlmm.nn.base import Dense
from torch import nn
import dgl

class message_v(nn.Module):
    def __init__(self,
                n_basis,
                n_atom_basis,
                activation=None,
                cutoff_network=None,
                cross_product = False) -> None:
        
        super(message_v,self).__init__()
        self.filter_vv = MLP(n_basis,n_atom_basis,n_hidden=n_basis,activation=activation)
        self.filter_vs = MLP(n_basis,n_atom_basis,n_hidden=n_atom_basis,activation=activation)
        self.transform_vv = MLP(n_atom_basis,n_atom_basis,n_hidden=n_atom_basis,activation=activation) 
        self.transform_vs = MLP(n_atom_basis,n_atom_basis,n_hidden=n_atom_basis,activation=activation) 
        self.cutoff_network = cutoff_network
        self.cross_product = cross_product

        if cross_product:
            self.cp1 = MLP(n_basis,n_atom_basis,n_hidden=n_atom_basis,activation=activation)
            self.cp2 = Dense(3,1,bias=False)
        self.coef= torch.nn.Parameter(torch.ones(n_atom_basis))
        self.bias= torch.nn.Parameter(torch.ones(n_atom_basis))
    
    def forward(self, g, feat_s, fij, vij, rij, feat_v=None):
        feat_v_norm = torch.norm(feat_v+1e-8,dim=-1)
        feat_v = feat_v/(((self.coef[None]**2)*feat_v_norm+self.bias[None])[:,:,None]+1e-7)

        if feat_v is not None:
            s_vv = self.transform_vv(feat_s)
            w_vv = self.filter_vv(fij)
            if self.cutoff_network is not None:
                C = self.cutoff_network(rij)
                w_vv = w_vv * C[:,None]
            m_vv = dgl.ops.e_mul_u(g,w_vv[...,None],s_vv[...,None]*feat_v)
            m_vv = dgl.ops.copy_e_sum(g,m_vv)
            
        s_vs = self.transform_vs(feat_s)

        if self.cross_product:
            w_vs1 = self.filter_vs(fij)
            w_vs2 = self.cp1(fij)
            if self.cutoff_network is not None:
                w_vs1 = w_vs1 * C[:,None]
                w_vs2 = w_vs2 * C[:,None]
            m_vs1 = dgl.ops.e_mul_u(g,vij[:,None]*w_vs1[...,None],s_vs[...,None])
            m_vs1 = dgl.ops.copy_e_sum(g,m_vs1)
            m_vs2 = dgl.ops.e_mul_u(g,vij[:,None]*w_vs2[...,None],s_vs[...,None])
            m_vs2 = dgl.ops.copy_e_sum(g,m_vs2)
            m_vs3 = torch.cross(m_vs1,m_vs2)
            m_vs = self.cp2(torch.stack([m_vs1,m_vs2,m_vs3],dim=-1)).squeeze()
        else:
            w_vs = self.filter_vs(fij)
            if self.cutoff_network is not None:
                w_vs = w_vs * C[:,None]
            m_vs = dgl.ops.e_mul_u(g,vij[:,None]*w_vs[...,None],s_vs[...,None])
            m_vs = dgl.ops.copy_e_sum(g,m_vs)

        m_v = m_vs + m_vv
        return m_v

class dipole_layer_eq(nn.Module):
    def __init__(self,n_atom_basis,activation=None):
        super(dipole_layer_eq,self).__init__()
        self.V = Dense(n_atom_basis,n_atom_basis,bias=False)
        self.mvv = MLP(2*n_atom_basis,n_atom_basis,activation=activation)

    def forward(self,feat_s, feat_v):
        feat_vv = self.V(feat_v.permute(0,2,1)).permute(0,2,1)
        feat_vv_norm = torch.norm(feat_vv+1e-9,dim=-1)
        input = torch.cat([feat_s,feat_vv_norm],dim=-1)
        mu_vv = self.mvv(input)
        #mu_vv = mu_vv[...,None]*(feat_vv/(feat_vv_norm[...,None]+1e-8)) ## new
        mu_vv = mu_vv[...,None]*feat_vv
        return mu_vv

class gated_update_v(nn.Module):
    def __init__(self,n_atom_basis,activation=None):
        super(gated_update_v,self).__init__()
        self.mss = MLP(2*n_atom_basis,n_atom_basis,activation=activation)
        self.msv = MLP(2*n_atom_basis,n_atom_basis,activation=activation)
        self.mvv = MLP(2*n_atom_basis,n_atom_basis,activation=activation)
        self.V = Dense(n_atom_basis,n_atom_basis,bias=False)
        self.U = Dense(n_atom_basis,n_atom_basis,bias=False)
        self.coef_vv = torch.nn.Parameter(torch.ones(n_atom_basis))
        self.coef_vu = torch.nn.Parameter(torch.ones(n_atom_basis))

    def forward(self, feat_s, feat_v):
        feat_vv = self.V(feat_v.permute(0,2,1)).permute(0,2,1)
        feat_vu = self.U(feat_v.permute(0,2,1)).permute(0,2,1)
        feat_vv_norm = torch.norm(feat_vv+1e-8,dim=-1)
        feat_vu_norm = torch.norm(feat_vu+1e-8,dim=-1)
        feat_vv = feat_vv/(((self.coef_vv[None]**2)*feat_vv_norm+1)[:,:,None]+1e-7)
        feat_vu = feat_vu/(((self.coef_vu[None]**2)*feat_vu_norm+1)[:,:,None]+1e-7)

        input = torch.cat([feat_s,feat_vv_norm],dim=-1)
        mu_ss = self.mss(input)
        mu_sv = self.msv(input)
        mu_vv = self.mvv(input)
        mu_sv_dot = torch.sum(feat_vv * feat_vu,dim=-1)
        mu_sv = mu_sv * mu_sv_dot
        mu_vv = mu_vv[...,None]*feat_vu

        mu_s = mu_sv + mu_ss
        return mu_s, mu_vv
