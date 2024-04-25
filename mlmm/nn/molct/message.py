import torch
import logging
from torch import nn
import dgl
from dgl import ops
from dgl import function as fn
from torch.cuda.amp import autocast
#from dgl.nn import EdgeWeightNorm
import numpy as np
from mlmm.nn.base import Dense
from mlmm.nn.molct.edge_embedding import Positional_Embedding, distance
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

__all__ = ["MPNN","EgoAttention","Transition"]


class MPNN(nn.Module):
    r"""Continuous-filter convolution block used in SchNet module.

    Args:
        n_in (int): number of input (i.e. atomic embedding) dimensions.
        ###n_filters (int): number of filter dimensions.
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
        n_atom_basis,
        n_gaussians,
        n_heads,
        filter_network,
        attention_network,
        cutoff,
        cutoff_network=None,
        activation=None,
        trainable_gaussians=False
    ):
        super(MPNN, self).__init__()
        
        ###self.filter_network = filter_network
        self.filter_network_linear, self.filter_network_res = filter_network
        
        # Perform Ego-Attention:
        self.attention_network = EgoAttention(n_atom_basis, n_heads=n_heads )
        self.cutoff_network = cutoff_network

        ### layer for compute distance
        self.distances = distance()

        # layer for expanding interatomic distances in a basis
        self.positional_embedding = Positional_Embedding(
            n_gaussians=n_gaussians,
            trainable_gaussians=trainable_gaussians,
            activation=activation,
            cutoff=cutoff,
            cutoff_network=cutoff_network,
            distance_expansion=None,                  
        )        

    @autocast(enabled=use_fp16)
    def forward(self, g, t, cell=None):
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
        logger = logging.getLogger('rl.'+__name__)

        self.g = g.local_var()
        r_ij, v_ij = self.distances(g,cell)
        logger.debug(f"r_ij:{r_ij.sum()}")

        ### Positional embedding:
        _, f_ij, q_id = self.positional_embedding(r_ij)
        logger.debug(f"f_ij:{f_ij.sum()}")

        # pass expanded interactomic distances through filter block
        W_in = 2. * f_ij - 1. ### rescale to [-1,1]
        ###W = self.filter_network(W) ### [N_b,N_a,1+N_n,n_atom_basis]; filter_net is a linear upscaling of f_ij   
        
        W = self.filter_network_linear(W_in) ### [N_b,N_a,1+N_n,n_atom_basis]; filter_net is a linear upscaling of f_ij
        if self.filter_network_res is not None:
            W = W + self.filter_network_res(W_in)
        
        # apply cutoff
        if self.cutoff_network is not None:
            C = self.cutoff_network(r_ij) ### [N_b,N_a,1+N_n]
            
        W = W*C[:,None] ### [N_b,N_a,1+N_n,1]
        
        # Perform MultiHead Attention:
        m = self.attention_network(self.g, t, W, q_id, cutoff_mask=C)

        logger.debug(f"m:{m.sum()}")
        return m, W 

__all__ = ["EgoAttention"]


class EgoAttention(nn.Module):
    r"""Ego Attention block used in 3DT module.

    Args:
        n_in (int): number of input (i.e. atomic embedding) dimensions.
        
    """

    def __init__(
        self,
        n_atom_basis,
        n_heads=8,        
    ):
        super(EgoAttention, self).__init__()
        
        # dense layer as mh_attention
        assert (n_atom_basis%n_heads==0), "Mismatch Head Numbers."
        n_per_head = n_atom_basis//n_heads
        
        self.n_heads = n_heads
        self.n_per_head = n_per_head
        
        self.mh_q = Dense(n_atom_basis, n_atom_basis, bias=False, activation=None)        
        self.mh_k = Dense(n_atom_basis, n_atom_basis, bias=False, activation=None)
        self.mh_v = Dense(n_atom_basis, n_atom_basis, bias=False, activation=None)
        self.mh_o = Dense(n_atom_basis, n_atom_basis, bias=False, activation=None) 
        
        self.layer_norm_in = nn.LayerNorm([n_atom_basis]) ###(input.size()[-1])

    @autocast(enabled=use_fp16)
    def forward(self, g, t, W, q_id, cutoff_mask=None):
        """Compute convolution block.

        Args:
            e (torch.Tensor): Element Embedding.
                with (N_b, N_a, n_atom_basis) shape.
            q (torch.Tensor): input representation/embedding of atomic environments
                with (N_b, N_a, n_atom_basis) shape.
            k (torch.Tensor): input representation/embedding of atomic environments
                with (N_b, N_a, N_n, n_atom_basis) shape.
            v (torch.Tensor): input representation/embedding of atomic environments
                with (N_b, N_a, N_n, n_atom_basis) shape.
            pairwise_mask: (N_b, N_a, N_n)
            cutoff_mask : (N_b, N_a, N_n, 1)

        Returns:
            torch.Tensor: block output with (N_b, N_a, n_out) shape.

        """

        logger = logging.getLogger('rl.'+__name__)
        
        n_heads = self.n_heads
        n_per_head = self.n_per_head
        
        ### ------------------###
        ### Ego Attention Transformer
        ### ------------------###
        logger.debug(f"W:{W.sum()}")
        logger.debug(f"t:{t.sum()}")
                
        with g.local_scope():
        ### Query:
            q = ops.e_mul_v(g, W, g.ndata['h1']) + t
            q = self.layer_norm_in(q)

            q_vec = self.mh_q(q[q_id])
            k_vec = self.mh_k(q)
            v_vec = self.mh_v(q) 

            ### Dot-product Attention:
            dot_product = (ops.u_mul_e(g,q_vec,k_vec).reshape(-1,n_heads,n_per_head)).sum(-1)

            if use_fp16:
                logit =  dot_product/ np.sqrt(self.n_per_head).astype(np.float16) 
            else:
                logit =  dot_product/ np.sqrt(self.n_per_head).astype(np.float32) 
            att_score  = ops.edge_softmax(g,logit,norm_by='src')

            ####
            att_score = att_score * cutoff_mask[:,None]

            with autocast(enabled=False):
                att_score = ops.e_div_u(g,att_score.float(),ops.copy_e_sum(g.reverse(), att_score))[:]

            g = g.reverse()
            g.edata['he'] = (v_vec.reshape([-1,n_heads,n_per_head]) * \
                            att_score[:,:,None]).reshape([-1,int(n_per_head*n_heads)])

            # Perform message-aggregation:
            g.update_all(fn.copy_e('he','m'),fn.sum('m','h2'))

            # Universal Affine Transform:
            m_agg = self.mh_o(g.ndata['h2'])
            
            logger.debug(f"m_agg:{m_agg.sum()}")


        return m_agg

class Transition(nn.Module):
    ### Added by Justin
    r"""Transition block for updating atomic embeddings.

    Args:
        n_atom_basis (int): number of features to describe atomic environments.

    """

    def __init__(
        self,
        n_atom_basis,
        activation=None,
    ):
        super(Transition, self).__init__()

        # filter block used in interaction block
        ###self.layer_norm_in = nn.LayerNorm([n_atom_basis]) ###(input.size()[-1])
        
        self.layer_norm_in = nn.Sequential(
            ###nn.Dropout(dropout_rate),
            nn.LayerNorm([n_atom_basis]), ###(input.size()[-1])
        )
        
        self.transition_network = nn.Sequential(
            Dense(n_atom_basis, n_atom_basis, activation=activation),
            Dense(n_atom_basis, n_atom_basis, activation=None),
        )
                                                       

    @autocast(enabled=use_fp16)
    def forward(self, x, v, t):
        """Compute interaction output.

        Args:
            x (torch.Tensor): input representation/embedding of atomic environments
                with (N_b, N_a, n_atom_basis) shape.
            r_ij (torch.Tensor): interatomic distances of (N_b, N_a, N_nbh) shape.
            neighbors (torch.Tensor): indices of neighbors of (N_b, N_a, N_nbh) shape.
            neighbor_mask (torch.Tensor): mask to filter out non-existing neighbors
                introduced via padding.
            f_ij (torch.Tensor, optional): expanded interatomic distances in a basis.
                If None, r_ij.unsqueeze(-1) is used.

        Returns:
            torch.Tensor: block output with (N_b, N_a, n_atom_basis) shape.

        """
        ### Time-Position Embedding:
        
        x = x + v
        x = self.layer_norm_in(x)
        
        x_t = self.transition_network(x)
        
        x = x + x_t
        ###x = self.layer_norm_out(x)

        return x
    
