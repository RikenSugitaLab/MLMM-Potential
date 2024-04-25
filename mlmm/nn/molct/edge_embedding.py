import torch
from torch import nn
import dgl
import dgl.ops as F
from torch.cuda.amp import autocast
from mlmm.nn.cutoff import CosineCutoff
import numpy as np
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

class distance(nn.Module):
    def __init__(self):
        super(distance,self).__init__()
 
    @autocast(enabled=use_fp16)
    def forward(self,g,cell=None):
 
        dis_vec = F.u_sub_v(g,g.ndata['xyz'],g.ndata['xyz'])
 
        if cell is not None:
           cell = dgl.broadcast_edges(g,cell)
           dis_vec[:,0] = torch.fmod(dis_vec[:,0],cell[:,0])
           dis_vec[:,1] = torch.fmod(dis_vec[:,1],cell[:,1])
           dis_vec[:,2] = torch.fmod(dis_vec[:,2],cell[:,2])
 
        dis = (torch.zeros([dis_vec.shape[0]]).to(dis_vec.device))

        dis[abs(dis_vec).sum(-1)!=0] = torch.norm(dis_vec[abs(dis_vec).sum(-1)!=0],dim=-1)
        return dis, dis_vec

class Positional_Embedding(nn.Module):
    ### Added by Justin
    r"""Relative positional embeddings.

    Args:
        n_hidden (int): number of hidden units in the FFMLP. Usually larger than n_atom_basis (recommend: 4*n_atom_basis).
        normalize_filter (bool, optional): if True, divide aggregated filter by number
            of neighbors over which convolution is applied.

    """

    def __init__(
        self,       
        n_gaussians,
        trainable_gaussians=False,
        cutoff=8.0,
        cutoff_network=CosineCutoff,
        distance_expansion=None,
        activation=None,        
    ):
        super(Positional_Embedding, self).__init__()
        
        # layer for expanding interatomic distances in a basis
        if distance_expansion is None:            
            ### Added by Justin:
            self.distance_expansion = LogNormalSmearing(
                np.log(0.1), np.log(cutoff), n_gaussians, trainable=trainable_gaussians
            )
        else:
            self.distance_expansion = distance_expansion
            
    @autocast(enabled=use_fp16)
    def forward(self, r_ij):
        """Compute interaction output.

        Args:
            r_ij (torch.Tensor): distance between any two vertices connected with edge 
                with (N_edge) shape.

        Returns:
            f_iij (torch.Tensor): edge representation by expansion with log gaussian.
            (N_edge,N_gaussians)

        """
        ### Position Embedding:
        
        ### Generate fij:
        q_id = (r_ij==0)
        r_ij[r_ij==0] = r_ij[r_ij==0] + 0.01 ### In case of overflow
        f_ij = self.distance_expansion(r_ij) ### Distance Embedding only
        return r_ij, f_ij, q_id        

class LogNormalSmearing(nn.Module):
    r"""Smear layer using a set of Log-Normal functions.

    Args:
        start (float, optional): log(minimal distance), :math:`\mu_0`.
        stop (float, optional): log(maximal distance), could be equal to log(cutoff), :math:`\mu_{N_g}`
        n_gaussians (int, optional): total number of Gaussian functions, :math:`N_g`.
        centered (bool, optional): If True, Gaussians are centered at the origin and
            the offsets are used to as their widths (used e.g. for angular functions).
        trainable (bool, optional): If True, widths and offset of Gaussian functions
            are adjusted during training process.

    """

    def __init__(
        self, start=np.log(0.1), stop=np.log(10.0), n_gaussians=64, centered=False, trainable=False,
        fix_start=np.log(0.1), fix_stop=np.log(8.0), fix_n_gaussians=32,
    ):
        super(LogNormalSmearing, self).__init__()
        # compute offset and width of Gaussian functions
            
        offset = torch.linspace(start, stop, n_gaussians)
        ###widths = torch.FloatTensor((offset[1] - offset[0]) * torch.ones_like(offset))
        
        fix_offset = torch.linspace(fix_start, fix_stop, fix_n_gaussians)
        widths = torch.FloatTensor((fix_offset[1] - fix_offset[0]) * torch.ones_like(offset))

        if trainable:
            self.width = nn.Parameter(widths)
            self.offsets = nn.Parameter(offset)
        else:
            self.register_buffer("width", widths)
            self.register_buffer("offsets", offset)
        self.centered = centered

    @autocast(enabled=use_fp16)
    def forward(self, distances):
        """Compute smeared-gaussian distance values.

        Args:
            distances (torch.Tensor): interatomic distance values of
                (N_b x N_at x N_nbh) shape.

        Returns:
            torch.Tensor: layer output of (N_b x N_at x N_nbh x N_g) shape.

        """
        if use_fp16:
            self.width = self.width.half()
            self.offsets = self.offsets.half()

        with autocast(enabled=False):
            return log_normal_smearing(
                distances, self.offsets, self.width, centered=self.centered
             )
### End of Justin's comment
    
@autocast(enabled=use_fp16)        
def log_normal_smearing(distances, offset, widths, centered=False):
    r"""Smear interatomic distance values using Gaussian functions.

    Args:
        distances (torch.Tensor): interatomic distances of (N_b x N_at x N_nbh) shape.
        offset (torch.Tensor): offsets values of Gaussian functions.
        widths: width values of Gaussian functions.
        centered (bool, optional): If True, Gaussians are centered at the origin and
            the offsets are used to as their widths (used e.g. for angular functions).

    Returns:
        torch.Tensor: smeared distances (N_b x N_at x N_nbh x N_g).

    """
    if not centered:
        # compute width of Gaussian functions (using an overlap of 1 STDDEV)
        with autocast(enabled=False):
            coeff = -0.5 / torch.pow(widths, 2)
        # Use advanced indexing to compute the individual components
            diff = torch.log( distances[:, None] ) - offset[None, :]
    else:
        # if Gaussian functions are centered, use offsets to compute widths
        coeff = -0.5 / torch.pow(widths, 2) 
        # if Gaussian functions are centered, no offset is subtracted
        diff = torch.log( distances[:, None] )
    # compute smear distance values

    with autocast(enabled=False):
        gauss = torch.exp(coeff * torch.pow(diff, 2))
    return gauss
