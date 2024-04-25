import torch
from torch import nn

from mlmm.nn import Dense


__all__ = ["MPNN"]


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
        filter_network,
        attention_network,
        cutoff_network=None,
        activation=None,
        
        relation_network=None,
        edge_vec=None,
    ):
        super(MPNN, self).__init__()
        
        self.filter_network = filter_network
        self.attention_network = attention_network
        self.cutoff_network = cutoff_network
        
        self.relation_network = relation_network
        self.edge_vec = edge_vec
        ### edge: [natom,natom,d]


    def forward(self, e, x, t, r_ij, neighbors, pairwise_mask, f_ij=None ):
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

        if f_ij is None:
            f_ij = r_ij.unsqueeze(-1)

        # pass expanded interactomic distances through filter block
        W = 2. * f_ij - 1. ### rescale to [-1,1]
        W = self.filter_network(W) ### [N_b,N_a,1+N_n,n_atom_basis]; filter_net is a linear upscaling of f_ij   
        
        if self.edge_vec is not None:
            ### edge_vec is a one-hot vector
            R = self.relation_network(self.edge_vec).unsqueeze(0) ### [1,natom,natom,nbasis]
            W = W+R
        
        
        # apply cutoff
        if self.cutoff_network is not None:
            C = self.cutoff_network(r_ij) ### [N_b,N_a,1+N_n]
            C = C.unsqueeze(-1) ### [N_b,N_a,1+N_n,1]
            
        W = W*C ### [N_b,N_a,1+N_n,1]
        
        # pass initial embeddings through Dense layer
        y = x ### [N_b,N_a,n_atom_basis]
        
        # reshape y for element-wise multiplication by W
        nbh_size = neighbors.size()
        nbh = neighbors.view(-1, nbh_size[1] * nbh_size[2], 1)
        nbh = nbh.expand(-1, -1, y.size(2))
        y = torch.gather(y, 1, nbh)
        y = y.view(nbh_size[0], nbh_size[1], nbh_size[2], -1) ### [N_b,N_a,N_n,n_atom_basis]
        
        # expand y
        y_ii = x.unsqueeze(-2)
        y = torch.cat( [y_ii,y], dim=-2 ) ### [N_b,N_a,1+N_n,n_atom_basis]
            
            
        # Perform MultiHead Attention:
        m = self.attention_network(e, x, y, t, W, pairwise_mask, cutoff_mask=C)
        # return: message : [N_b,N_a,n_atom_basis]

        
        return m, W
