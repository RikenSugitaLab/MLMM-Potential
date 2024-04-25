import torch
from torch import nn
import dgl
from torch.nn.modules import activation
from mlmm.nn.base import Dense
from mlmm.nn.initializers import zeros_initializer
from mlmm.nn.activations import swish, swish2, shifted_softplus
from mlmm.nn.cutoff import *
import random 
import numpy as np
from typing import Union
from mlmm.nn import aggregator

__all__=["Dipole_Layer","FieldInteraction","TensorInteraction", "ChargeFieldInteraction","gvp"]

class Dipole_Layer(nn.Module):
    """
       Args:
        atom_features: 
        dipole_features:
        transform:
        cutoff:
        activation:
    """
    def __init__(self,
                atom_features: int,
                dipole_features: int,
                transform: bool = True,
                cutoff_network = None,
                activation = swish,
                cross_product = False
                ):
        super(Dipole_Layer, self).__init__()
        self.cross_product = cross_product

        if transform:
            self.transform = nn.Sequential(
                Dense(atom_features, atom_features, activation=activation),
                Dense(atom_features, dipole_features, activation=activation)
            )

            if cross_product:
                self.transform2 = nn.Sequential(
                Dense(atom_features, atom_features, activation=activation),
                Dense(atom_features, dipole_features, activation=activation)
            )
                self.mix = Dense(3,1, activation=None)
        else:
            self.transform = None
        self.cutoff_network = cutoff_network

    def forward(self, g, x, rij, vij):
        
        # Apply transformation layer
        if self.transform is not None:
            q = self.transform(x)
            if self.cross_product:
                q2 = self.transform2(x)
        else:
            q = x

        mu_ij = dgl.ops.e_mul_u(g, vij[:,None], q[:,:,None])


        if self.cutoff_network is not None:
            c_ij = self.cutoff_network(rij)
            mu_ij = mu_ij * c_ij[:, None, None]

            if self.cross_product:
                mu_ij2 = dgl.ops.e_mul_u(g, vij[:,None], q2[:,:,None])
                mu_ij2 = mu_ij2 * c_ij[:, None, None]

        # Form final sum
        mu = dgl.ops.copy_e_sum(g, mu_ij)
        if self.cross_product:
            mu2 = dgl.ops.copy_e_sum(g, mu_ij2)
            mu3 = torch.cross(mu, mu2)
            mu = self.mix(torch.stack([mu,mu2,mu3],dim=-1)).squeeze(-1)

        return mu

class TensorInteraction(nn.Module):

    def __init__(
            self,
            dipole_features: int,
            n_atom_basis: int,
            cutoff: float = 5.0,
            n_basis: int = 50,
            shielding=None,
            cutoff_function=MollifierCutoff,
            activation = swish,
            quadrupole = False,
            aggregator = aggregator):
        super(TensorInteraction, self).__init__()

        self.cutoff = cutoff_function(cutoff)
        self.shielding = shielding
        self.quadrupole = quadrupole
        self.aggregator = aggregator

        self.denseq = nn.Sequential(
            Dense(dipole_features, dipole_features, activation=activation),
            Dense(dipole_features, n_atom_basis)
        )

        self.denseq = nn.Sequential(
            Dense(n_atom_basis, n_atom_basis, activation=activation),
            Dense(n_atom_basis, n_atom_basis)
        )

        self.dense_v = Dense(dipole_features,n_atom_basis,bias=False)

        self.dense = nn.Sequential(
            Dense(dipole_features, dipole_features, activation=activation),
            Dense(dipole_features, n_atom_basis)
        )

        self.distance_expansion = nn.Sequential(
            Dense(n_basis, dipole_features, activation=activation),
            Dense(dipole_features, dipole_features, bias_init=zeros_initializer)
        )

    def forward(self, g, mu, rij, vij, f_ij, feat):

        # Matrix multiplication with interaction tensor expressed as dot products
        # mu_i @ Tij @ mu_j.T = mu_i @ ( -1 rij^2 + 3 Rij.T Rij ) @ mu_j
        #                     = -(mu_i @ mu_j.T)*rij^2 + 3 (mu_i @ Rij.T)*(mu_j @ Rij)
        diagonal_term =  torch.sum(dgl.ops.u_mul_v(g, mu, mu), dim=-1) * (rij**2)[:,None] ## need test
        outer_term1 = torch.sum(dgl.ops.u_mul_e(g, mu, vij[:,None]), dim=-1)
        outer_term2 = torch.sum(dgl.ops.v_mul_e(g, mu, vij[:,None]), dim=-1)
        outer_term = outer_term1 * outer_term2

        # Inverse r^5 in this term
        radial = self.distance_expansion(f_ij) / ((rij[..., None] ** 5)+1e-10)

        # Apply shielding
        if self.shielding is not None:
            #quadrupole = dgl.ops.copy_e_sum(g,quadrupole)
            radial = radial * (1 - self.shielding(rij)[..., None])

        if self.cutoff is not None:
            radial = radial * self.cutoff(rij)[..., None]

        if "T2" in g.ndata.keys() and "T3" in g.ndata.keys():
            q = self.denseq(feat)
            outer_T = vij[:,:,None]@vij[:,None]
            diag_T= torch.eye(3)[None].to(rij.device)*(rij**2)[:,None,None]
            T = diag_T - 3*outer_T
            T = T[:,None] * radial[:,:,None,None]
            mv1 = dgl.ops.e_mul_u(g,dgl.ops.e_dot_v(g,T,mu[:,:,None]).squeeze(),q[...,None])
            mv2 = dgl.ops.e_mul_v(g,dgl.ops.e_dot_u(g,T,mu[:,:,None]).squeeze(),q[...,None])
            mvt = dgl.ops.copy_e_sum(g, mv1-mv2)
            mvt = self.dense_v(mvt.permute(0,2,1)).permute(0,2,1)
        else:
            mvt = 0

        v = (diagonal_term - 3 * outer_term) * radial

        # Sum over neighbors
        # v = dgl.ops.copy_e_sum(g, v)
        v = self.aggregator(g, v)

        # Apply final transformation
        v = self.dense(v)

        return v, mvt


class multipole_interaction(nn.Module):
    def __init__(self,
                 dipole_features: int, 
                 n_atom_basis: int, 
                 activation = swish):
        super(multipole_interaction,self).__init__()

        self.denseq = nn.Sequential(
                        Dense(dipole_features, dipole_features, activation=activation),
                        Dense(dipole_features, dipole_features))
        self.dense_T1 = Dense(dipole_features,dipole_features,bias=False)
        self.dense_T2 = Dense(dipole_features,dipole_features,bias=False)
        self.dense_T3 = Dense(dipole_features,dipole_features,bias=False)

        self.dense_s = nn.Sequential(
                        Dense(dipole_features, dipole_features, activation=activation),
                        Dense(dipole_features, n_atom_basis))
        self.dense_v = Dense(dipole_features,n_atom_basis,bias=False)
    
    def forward(self, g, feat, mu):
        M0 = self.denseq(feat)

        M1 = self.dense_T1(mu.permute(0,2,1)).permute(0,2,1)
        ms = torch.sum(M1*g.ndata['T1'].view(-1,1,3),dim=-1)
        ms = ms + (g.ndata['T0'][:,None]*M0)

        mv = g.ndata['T1'][:,None]*M0[...,None]

        if 'T2' in g.ndata.keys() and 'T3' in g.ndata.keys():
            M2 = self.dense_T2(mu.permute(0,2,1)).permute(0,2,1)
            quadrupole = M2[...,None]@M2[:,:,None]
            ms = ms + (quadrupole * g.ndata['T2'][:,None]).sum(-1).sum(-1)
            mv = mv + torch.matmul(g.ndata['T2'][:,None],M1[...,None]).squeeze()

        if 'T3' in g.ndata.keys() and 'T4' in g.ndata.keys():
            M3 = self.dense_T3(mu.permute(0,2,1)).permute(0,2,1)
            octupole = M3[:,:,:,None,None]*M3[:,:,None,:,None]*M3[:,:,None,None,:]
            ms = ms + (octupole * g.ndata['T3'][:,None]).sum(-1).sum(-1).sum(-1)
            mv = mv + (g.ndata['T3'][:,None]*quadrupole[:,:,None]).sum(-1).sum(-1)

        mv = self.dense_v(mv.permute(0,2,1)).permute(0,2,1)
        ms = self.dense_s(ms)

        return ms, mv




        
