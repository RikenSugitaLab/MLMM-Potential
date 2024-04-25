import torch
from torch import nn
import dgl
import logging
import dgl.ops as F
import dgl.function as fn
from mlmm.nn.activations import shifted_softplus, swish
from mlmm.nn.aggregator import aggregator
from mlmm.utils.utils import check_value
from mlmm.nn.base import Dense
from mlmm.nn.attension import attension
from mlmm.nn.acsf import GaussianSmearing, LogGaussianDistribution 
from torch.cuda.amp import autocast
import os

__all__=["distance","compute_message"]
class distance(nn.Module):
    def __init__(self):
        super(distance,self).__init__()
    
    def forward(self,g,cell=None):
        logger = logging.getLogger('rl.'+__name__)

        dis_vec = F.u_sub_v(g,g.ndata['xyz'],g.ndata['xyz'])

        if cell is not None:
           cell = dgl.broadcast_edges(g,cell)
           #dis_vec[:,0] = torch.fmod(dis_vec[:,0],cell[:,0])
           #dis_vec[:,1] = torch.fmod(dis_vec[:,1],cell[:,1])
           #dis_vec[:,2] = torch.fmod(dis_vec[:,2],cell[:,2])
           dis_vec2 = torch.cat([(cell[:,:3] - abs(dis_vec))[:,:,None],abs(dis_vec)[:,:,None]],dim=-1)
           dis_vec2 = dis_vec2.min(-1)[0]
        mask2 = torch.where(abs(dis_vec2)==abs(dis_vec),1,-1)
        mask = (dis_vec+1e-8)/abs(dis_vec+1e-8)
        mask = mask2 * mask

        dis = torch.norm(dis_vec+1e-9,dim=-1)
        if check_value(dis_vec):
            logger.critical('value error: nan or inf')
            assert check_value(dis_vec)==0 

        if check_value(mask):
            logger.critical('value error: nan or inf')
            assert check_value(mask)==0 
        return dis, dis_vec2 * mask 
          
class compute_message(nn.Module):
    r"""Continuous-filter convolution block used in SchNet module.

    Args:
        n_in (int): number of input (i.e. atomic embedding) dimensions.
        n_filters (int): number of filter dimensions.
        n_out (int): number of output dimensions.
        filter_network (nn.Module): filter block.
        cutoff_network (nn.Module, optional): if None, no cut off function is used.
        activation (callable, optional): if None, no activation function is used.
        normalize_filter (bool, optional): If True, normalize filter to the number
            of neighbors when aggregating.
        axis (int, optional): axis over which convolution should be applied.

    """

    def __init__(
        self,
        n_in,
        n_filters,
        filter_network,
        cutoff_network=None,
        activation=swish,
        aggregator=aggregator
    ):
        super(compute_message, self).__init__()
        self.in2f = Dense(n_in, n_filters, bias=False, activation=None)
        self.filter_network = filter_network
        self.cutoff_network = cutoff_network

        # self.att = attension(n_in,n_filters)
        self.msg_transform = nn.Sequential(Dense(n_filters, n_in, bias=True, activation=activation),
                             Dense(n_in,n_in,bias=True))
        self.aggregator = aggregator

    def forward(self, g, feat, fij, rij): 
        """Compute convolution block.

        Args:
            x (torch.Tensor): input representation/embedding of atomic environments
                with (N_b, N_a, n_in) shape.
            r_ij (torch.Tensor): interatomic distances of (N_b, N_a, N_nbh) shape.
            neighbors (torch.Tensor): indices of neighbors of (N_b, N_a, N_nbh) shape.
            pairwise_mask (torch.Tensor): mask to filter out non-existing neighbors
                introduced via padding.
            f_ij (torch.Tensor, optional): expanded interatomic distances in a basis.
                If None, r_ij.unsqueeze(-1) is used.

        Returns:
            torch.Tensor: block output with (N_b, N_a, n_out) shape.

        """
        ### continuous filter ###
        with g.local_scope():

            logger = logging.getLogger('rl.'+__name__)
            logger.debug(f'fij: {fij.sum()}')
            if check_value(fij):
                logger.critical('value error: nan or inf')
                assert check_value(fij)==0 

        # apply cutoff
            if self.cutoff_network is not None:
                C = self.cutoff_network(rij)
                logger.debug(f"C:{C.sum()}")
                g.edata['he'] = self.filter_network(fij) * C.unsqueeze(-1)#*self.att(g,feat)
            else:
                g.edata['he'] = self.filter_network(fij) #* self.att(g,feat)

            logger.debug(f"feat:{feat.sum()}")
            g.ndata['hv'] = self.in2f(feat)
            m = dgl.ops.u_mul_e(g,g.ndata['hv'],g.edata['he'])
            m = self.aggregator(g,m)

        ### compute message ###

            return self.msg_transform(m.squeeze())
