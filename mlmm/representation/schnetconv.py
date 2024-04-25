import dgl
import torch as th
import dgl.ops as F
import dgl.function as fn
from mlmm.nn.aggregator import aggregator
from mlmm.nn.cfconv import compute_message, distance
from mlmm.nn.cutoff import CosineCutoff, MollifierCutoff
import logging
from mlmm.nn.base import Dense
from mlmm.nn.activations import swish,swish2,shifted_softplus
from torch import nn
from torch.autograd import grad
from mlmm.nn.graph_norm import graph_norm,gen_read_out
from mlmm.nn.blocks import MLP
from dgl.nn.pytorch.glob import SumPooling
from mlmm.utils.logger import setup_logger
from mlmm.utils.utils import check_value,add_sn
from mlmm.nn.acsf import GaussianSmearing, LogGaussianDistribution 
import os

__all__=["SchNetConv,SchNet"]

class SchNetConv(nn.Module):
    def __init__(self,
                 n_atom_basis=64,
                 n_filters=64,
                 n_basis=25,
                 cutoff_network=CosineCutoff,
                 activation=swish,
                 aggregator = aggregator,
                 ):
        
        super(SchNetConv,self).__init__()
        self.g = None
        self.filter_network = nn.Sequential(
            Dense(n_basis, n_filters, activation=activation),
            Dense(n_filters, n_filters),
        )

        # cutoff layer used in interaction block
        self.cutoff_network = cutoff_network
        self.message = compute_message(n_atom_basis,
                                       n_filters,
                                       self.filter_network,
                                       cutoff_network=self.cutoff_network,
                                       activation=activation,
                                       aggregator=aggregator
                                      )
        self.n_atom_basis, self.n_filters, self.n_basis= (n_atom_basis,n_filters,n_basis)
        self.activation = activation

    def forward(self, g, feat , fij, rij):
        logger = logging.getLogger('rl.'+__name__)
        logger.debug(f"schnetconv start")
        self.g = g.local_var()

        m = self.message(self.g, feat ,fij,rij)

        logger.debug(f"message: {m.sum()}")
        if check_value(m):
            logger.critical('value error: nan or inf')
            assert check_value(m)==0 

        return m 

class SchNet(nn.Module):
    def __init__(self,
                 n_interactions=3,
                 n_atom_basis=64,
                 n_filters=64,
                 n_basis=25,
                 cutoff=5,
                 max_z=100,
                 cutoff_network: str = 'CosineCutoff',
                 activation: str = 'swish',
                 graph_rep = False,
                 basis_func = 'gaussian',
                 ):
        super(SchNet,self).__init__()

        self.n_atom_basis = n_atom_basis
        self.embedding = nn.Embedding(max_z, n_atom_basis, padding_idx=0)
        self.layers = nn.ModuleList()
        self.graph_norm = nn.ModuleList() 

        self.distance = distance()

        if cutoff_network == "MollifierCutoff":
            cutoff_network = MollifierCutoff(cutoff=cutoff)
        elif cutoff_network == "CosineCutoff":
            cutoff_network = CosineCutoff(cutoff=cutoff)

        if basis_func == 'loggaussian':
             self.distance_expansion = LogGaussianDistribution(d_max=cutoff,num_rbf=n_basis)
        elif basis_func == 'gaussian':
             self.distance_expansion = GaussianSmearing(
                 0.0, cutoff, n_gaussians=n_basis)

        # if spectral_norm:
            # self.out_net = nn.Sequential(
                    # MLP(n_atom_basis, 1, None, 2, swish2),
                # )
            # self.out_net = add_sn(self.out_net)
        # else:
            # self.out_net = nn.Sequential(
                    # MLP(n_atom_basis, 1, None, 2, swish),
                # )


        self.graph_rep = graph_rep
        if graph_rep:
            self.rep_pooling = gen_read_out(agg='softmax')

        for _ in range(n_interactions):
            self.layers.append(SchNetConv(
                                          n_atom_basis=n_atom_basis,
                                          n_filters=n_filters,
                                          n_basis=n_basis,
                                          cutoff=cutoff,
                                          cutoff_network=cutoff_network))

            self.graph_norm.append(graph_norm(n_atom_basis))

    def forward(self,g, cell=None):
        logger = logging.getLogger('rl.'+__name__)
        self.g = g
        
        ## calculate distance and distance expansion
        rij, vij = self.distance(g,cell=cell)
        fij = self.distance_expansion(rij)

        if (rij>self.cutoff+0.1).sum()>0:
            logger.debug(f"{(rij>self.cutoff).sum()}")
            logger.debug(f"{rij[rij>self.cutoff]}")
            logger.critical('error: distance calculation errror, distance should not be large than cutoff distance')
            assert (rij>self.cutoff).sum()==0

        feat = self.embedding(self.g.ndata['z'].long())

        self.feat_list = []
        self.feat_list.append(feat)

        for layer, gn in zip(self.layers, self.graph_norm):
            feat = layer(self.g,feat,fij,rij,cell) + feat
            feat = gn(self.g, feat)
            self.feat_list.append(feat)

        return feat, self.rep_pooling(self.g,self.feat_list[-1]) 

        
        
        
        



        
