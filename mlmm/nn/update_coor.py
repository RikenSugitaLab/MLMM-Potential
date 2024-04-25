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
__all__=['Update_coor']

class Update_coor(nn.Module):
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
        super(Update_coor, self).__init__()
        self.cross_product = cross_product

        if transform:
            self.transform = nn.Sequential(
                Dense(atom_features, atom_features, activation=activation),
                Dense(atom_features, 1, activation=activation)
            )

            if cross_product:
                self.transform2 = nn.Sequential(
                Dense(atom_features, atom_features, activation=activation),
                Dense(atom_features, dipole_features, activation=activation)
            )
                self.mix = Dense(3,1, activation=None)
        else:
            self.transform = None
        self.comb = Dense(dipole_features,1,bias=False)
        self.cutoff_network = cutoff_network

    def forward(self, g, x, rij, vij):
        
        # Apply transformation layer
        if self.transform is not None:
            q = self.transform(x)
            if self.cross_product:
                q2 = self.transform2(x)
        else:
            q = x

        mu_ij = dgl.ops.e_mul_u(g, vij, q)


        if self.cutoff_network is not None:
            c_ij = self.cutoff_network(rij)
            mu_ij = mu_ij * c_ij[:, None]

            if self.cross_product:
                mu_ij2 = dgl.ops.e_mul_u(g, vij, q2)
                mu_ij2 = mu_ij2 * c_ij[:, None]
        
        # Form final sum
        mu = dgl.ops.copy_e_sum(g, mu_ij)
        if self.cross_product:
            mu2 = dgl.ops.copy_e_sum(g, mu_ij2)
            mu3 = torch.cross(mu, mu2)
            mu = self.mix(torch.stack([mu,mu2,mu3],dim=-1)).squeeze(-1)

        return mu.squeeze()
