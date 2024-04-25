import torch
from torch import nn
from torch.cuda.amp import autocast
from mlmm.nn.cutoff import CosineCutoff
from mlmm.nn.blocks import MLP
from mlmm.nn.activations import swish
import math
from scipy import special
import dgl
from mlmm.nn.base import Dense

__all__ = ["gaussian_smearing","GaussianSmearing", "LogGaussianDistribution", "GTO_Expansion" ]
def gaussian_smearing(distances,offset, widths, centered=False):
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
        coeff = -0.5 / torch.pow(widths, 2)
        # Use advanced indexing to compute the individual components
        diff = distances[:, None] - offset[None, :]
    else:
        # if Gaussian functions are centered, use offsets to compute widths
        coeff = -0.5 / torch.pow(offset, 2)
        # if Gaussian functions are centered, no offset is subtracted
        diff = distances[:, None]
    # compute smear distance values
    gauss = torch.exp(coeff * torch.pow(diff, 2))
    return gauss


class GaussianSmearing(nn.Module):
    r"""Smear layer using a set of Gaussian functions.

    Args:
        start (float, optional): center of first Gaussian function, :math:`\mu_0`.
        stop (float, optional): center of last Gaussian function, :math:`\mu_{N_g}`
        n_gaussians (int, optional): total number of Gaussian functions, :math:`N_g`.
        centered (bool, optional): If True, Gaussians are centered at the origin and
            the offsets are used to as their widths (used e.g. for angular functions).
        trainable (bool, optional): If True, widths and offset of Gaussian functions
            are adjusted during training process.

    """

    def __init__(
        self, start=0.0, stop=5.0, n_gaussians=50, centered=False, trainable=False
    ):
        super(GaussianSmearing, self).__init__()
        # compute offset and width of Gaussian functions
        offset = torch.linspace(start, stop, n_gaussians)
        # widths = torch.FloatTensor((offset[1] - offset[0]) * torch.ones_like(offset))
        widths = (offset[1] - offset[0]) * torch.ones_like(offset)

        if trainable:
            self.width = nn.Parameter(widths)
            self.offsets = nn.Parameter(offset)
        else:
            self.register_buffer("width", widths)
            self.register_buffer("offsets", offset)
        self.centered = centered

    def forward(self,rij):
        """Compute smeared-gaussian distance values.

        Args:
            distances (torch.Tensor): interatomic distance values of
                (N_b x N_at x N_nbh) shape.

        Returns:
            torch.Tensor: layer output of (N_b x N_at x N_nbh x N_g) shape.

        """
        if rij.dtype==torch.float16:
           self.width = self.width.half()
           self.offsets = self.offsets.half()

        with autocast(enabled=False):
            return gaussian_smearing(
                rij, self.offsets, self.width, centered=self.centered
            )

class LogGaussianDistribution(nn.Module):
    def __init__(
        self,
        d_min=1e-3,
        d_max=1.0,
        num_rbf=32,
        sigma=None,
        trainable=False,
        min_cutoff=False,
        max_cutoff=False,
    ):
        super().__init__()
        if d_max <= d_min:
            raise ValueError('The argument "d_max" must be larger'+
                'than the argument "d_min" in LogGaussianDistribution!')
            
        if d_min <= 0:
            raise ValueError('The argument "d_min" must be '+
                ' larger than 0 in LogGaussianDistribution!')
            
        self.d_max = d_max
        self.d_min = d_min / d_max
        self.min_cutoff=min_cutoff
        self.max_cutoff=max_cutoff
        
        log_dmin=math.log(self.d_min)
        centers = torch.linspace(log_dmin,0,num_rbf)
        self.ones = torch.ones_like(centers)

        if sigma is None:
            sigma = -log_dmin / (num_rbf-1)
        rescale = (-0.5 / (sigma * sigma))*self.ones

        if trainable:
            self.rescale = nn.Parameter(rescale)
            self.centers = nn.Parameter(centers)
        else:
            self.register_buffer("rescale", rescale)
            self.register_buffer("centers", centers)


    def forward(self,distance):

        dis = distance / self.d_max
       
        if self.min_cutoff:
            dis = torch.max(dis,self.d_min)

        log_dis = torch.log(dis[:,None])
        log_diff = log_dis - self.centers
        log_diff2 = torch.pow(log_diff,2)
        log_gauss = torch.exp( self.rescale * log_diff2  )

        if self.max_cutoff:
            ones = torch.ones_like(exdis)
            zeros = torch.zeros_like(exdis)
            cuts = F.select(exdis < 1.0, ones, zeros)
            log_gauss = log_gauss * cuts
    
        return log_gauss

class GTO_Expansion(torch.nn.Module):
    def __init__(self, 
                start: float,
                stop: float,
                n_basis: int,
                n_atom_basis: int,
                L:int,
                cutoff_network = None,
                activation = swish
                ):
        super(GTO_Expansion,self).__init__()
        offset = torch.linspace(start, stop, n_basis)
        width = (offset[1] - offset[0]) * torch.ones_like(offset)
        self.register_buffer("width", width)
        self.register_buffer("offsets", offset)
        self.cutoff_network = cutoff_network

        self.L = L
        self.i = []
        self.j = []
        self.k = []
        for i in range(L+1):
            for j in range(L+1-i):
                k = L - i - j
                self.i.append(i)
                self.j.append(j)
                self.k.append(k)
        self.i = torch.FloatTensor(self.i)
        self.j = torch.FloatTensor(self.j)
        self.k = torch.FloatTensor(self.k)
        self.register_buffer("ipower", self.i)
        self.register_buffer("jpower", self.j)
        self.register_buffer("kpower", self.k)

        d = special.factorial(self.i)*special.factorial(self.j)*special.factorial(self.k)
        prefactor2 = torch.sqrt(special.factorial(L)/d)
        self.register_buffer("prefactor2", prefactor2)

        self.cj_mlp = MLP(n_atom_basis, 1, activation=activation)
        self.msg = nn.Sequential(
            Dense(n_basis, n_atom_basis, activation=activation),
            Dense(n_atom_basis, n_atom_basis),
        )

        self.filter = nn.Sequential(
            Dense(self.i.shape[0], self.i.shape[0], activation=activation,bias=False),
            Dense(self.i.shape[0], 1),
        )

    def forward(self, g, feat, dis_vec):
        phix = torch.pow(dis_vec[:,0][:,None]+1e-8, self.ipower[None,:])
        phiy = torch.pow(dis_vec[:,1][:,None]+1e-8, self.jpower[None,:])
        phiz = torch.pow(dis_vec[:,2][:,None]+1e-8, self.kpower[None,:])
        prefactor = phix*phiy*phiz

        dis = torch.norm(dis_vec+1e-9, dim=-1)
        if self.cutoff_network is not None:
            C = self.cutoff_network(dis)
        else: 
            C = 1

        # compute width of Gaussian functions (using an overlap of 1 STDDEV)
        coeff = -0.5 / torch.pow(self.width, 2)
        # Use advanced indexing to compute the individual components
        diff = dis[:, None] - self.offsets[None, :]
        gauss = C[:,None]*torch.exp(coeff * torch.pow(diff, 2))

        gauss = prefactor[:,None]*gauss[...,None]*self.prefactor2[None,None]

        cj = self.cj_mlp(feat)
        fij = dgl.ops.e_mul_v(g, gauss, cj)
        fij = dgl.ops.copy_e_sum(g.reverse(), fij)
        fij = (fij**2).sum(-1)

        #msg= fij/(torch.norm(fij+1e-9,dim=-1)[:,None] + 1e-9)
        msg= fij/(torch.norm(fij+1e-9,dim=-1)[:,None] + 1.0)

        #return self.msg(msg),self.filter(gauss/(torch.norm(gauss+1e-8,dim=-1)+1e-8)[...,None]).squeeze()
        return self.msg(msg),self.filter(gauss/(torch.norm(gauss+1e-8,dim=-1)+1.0)[...,None]).squeeze()
