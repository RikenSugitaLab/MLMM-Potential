import torch
from torch import nn
from mlmm.nn.activations import swish
from mlmm.nn.cutoff import CosineCutoff
import dgl
from dgl.nn.pytorch.glob import SumPooling
from torch.cuda.amp import autocast
from mlmm.nn.molct.edge_embedding import distance, Positional_Embedding
from mlmm.nn.base import Dense
from mlmm.nn.molct.message import MPNN, EgoAttention, Transition
from mlmm.nn.act import AdaptiveComputationTime, _gen_timing_signal
import numpy as np
import logging
from mlmm.nn.blocks import MLP
import os

__all__ = ["TDT_Interaction","TDTNet"]

class TDT_Interaction(nn.Module):
    r"""SchNet interaction block for modeling interactions of atomistic systems.

    Args:
        n_atom_basis (int): number of features to describe atomic environments.
        n_spatial_basis (int): number of input features of filter-generating networks.
        n_gaussians (int): number of log-gaussians used to expand distance
        n_heads (int): number of heads in transformer
        cutoff (float): cutoff radius.
        cutoff_network (nn.Module, optional): cutoff layer.
        activation (nn.Module)
        resnet(bool): if use resnet
        post_ln(bool): if layer normalization
        apply_transition_function(bool): if use transition function

    """

    def __init__(
        self,
        n_atom_basis,
        n_gaussians,
        n_heads,

        cutoff=8.0,        
        cutoff_network=CosineCutoff,
        activation=None,

        apply_transition_function=False,  
        post_ln=False,
        
        resnet=True,
    ):
        super(TDT_Interaction, self).__init__()

        # filter block used in interaction block 
        ##self.filter_network = Dense(n_gaussians, n_atom_basis, bias=True, activation=None) 
        self.filter_network_linear = Dense(n_gaussians, n_atom_basis, bias=True, activation=None) 
        if resnet:
            self.filter_network_res = nn.Sequential(
                Dense(n_gaussians, n_atom_basis, activation=activation),
                Dense(n_atom_basis, n_atom_basis, activation=None),
            )
            self.filter_network = [self.filter_network_linear,self.filter_network_res]
        else:
            self.filter_network = [self.filter_network_linear,None]
    
        # cutoff layer used in interaction block
        self.cutoff_network = cutoff_network(cutoff)
        self.attention_network = EgoAttention(n_atom_basis, n_heads=n_heads )
        
        # For Message Passing (interaction block)
        self.mpnn = MPNN(
            n_atom_basis,
            n_gaussians,
            n_heads,            
            self.filter_network, ### Filter_Network is for positional embedding
            self.attention_network,
            cutoff,
            cutoff_network=self.cutoff_network,
            activation=activation,
        )
        
        # Transition function:
        self.apply_transition_function = apply_transition_function
        if apply_transition_function:
            self.transition = Transition(n_atom_basis=n_atom_basis, activation=activation)  
            
        self.post_ln = post_ln
        if post_ln:
            self.layer_norm = nn.LayerNorm([n_atom_basis])
            
    def forward(self, g, x, t, cell=None):
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
        # continuous-filter convolution interaction block followed by Dense layer
        self.g = g.local_var()
        self.g.ndata['h1'] = x
        v,W = self.mpnn(self.g, t, cell=cell)

        if self.apply_transition_function:
            x = self.transition(x, v, t, W) 
        else:
            x = x + v
            
        ###x = x + v
        if self.post_ln:
            x = self.layer_norm(x)
        
        return x

########################################################

class TDTNet(nn.Module):
    """SchNet architecture for learning representations of atomistic systems.

    Args:
        n_atom_basis (int, optional): number of features to describe atomic environments.
            This determines the size of each embedding vector; i.e. embeddings_dim.
        n_interactions (int, optional): number of interaction blocks.
        n_scales (int, optional): number of interaction blocks.
        cutoff (float, optional): cutoff radius.
        n_gaussians (int, optional): number of Gaussian functions used to expand
            atomic distances.
        normalize_filter (bool, optional): if True, divide aggregated filter by number
            of neighbors over which convolution is applied.
        resnet(bool): if use resnet
        post_ln(bool): if layer normalization
        apply_transition_function(bool): if use transition function
        use_act(bool): if use Adaptive Computation Time
        coupled_interactions (bool, optional): if True, share the weights across
            interaction blocks and filter-generating networks.
        return_intermediate (bool, optional): if True, `forward` method also returns
            intermediate atomic representations after each interaction block is applied.
        max_z (int, optional): maximum nuclear charge allowed in database. This
            determines the size of the dictionary of embedding; i.e. num_embeddings.
        cutoff_network (nn.Module, optional): cutoff layer.
        trainable_gaussians (bool, optional): If True, widths and offset of Gaussian
            functions are adjusted during training process.
        distance_expansion (nn.Module, optional): layer for expanding interatomic
            distances in a basis.
        charged_systems (bool, optional):

    References:
    .. [#schnet1] Schütt, Arbabzadah, Chmiela, Müller, Tkatchenko:
       Quantum-chemical insights from deep tensor neural networks.
       Nature Communications, 8, 13890. 2017.
    .. [#schnet_transfer] Schütt, Kindermans, Sauceda, Chmiela, Tkatchenko, Müller:
       SchNet: A continuous-filter convolutional neural network for modeling quantum
       interactions.
       In Advances in Neural Information Processing Systems, pp. 992-1002. 2017.
    .. [#schnet3] Schütt, Sauceda, Kindermans, Tkatchenko, Müller:
       SchNet - a deep learning architecture for molceules and materials.
       The Journal of Chemical Physics 148 (24), 241722. 2018.

    """

    def __init__(
        self,
        ### Newtork Hyper-parameters:
        n_atom_basis=128,
        n_gaussians=32, ### 25
        n_heads=8,
        
        activation=swish,
        
        ### Model Hyper-parameters:
        n_interactions=4,
        n_scales=1,        
        cutoff=8.0,
        
        apply_transition_function=False, ### If true, Apply Transition function as in Transformer
        post_ln=False,

        use_act=True, ### Adaptive Computation Time        
        
        return_intermediate=False,
        max_z=100,
        cutoff_network=CosineCutoff,
        
        distance_expansion=None,
        charged_systems=False,
        
        resnet = True,
    ):
        super(TDTNet, self).__init__()

        self.n_atom_basis = n_atom_basis
        # make a lookup table to store embeddings for each element (up to atomic
        # number max_z) each of which is a vector of size n_atom_basis
        self.embedding = nn.Embedding(max_z, n_atom_basis, padding_idx=0)

        # layer for computing interatomic distances
        self.distances = distance()

        
        # block for computing interaction
        
        self.interaction_blocks = nn.ModuleList(
                [
                    TDT_Interaction(
                        n_atom_basis=n_atom_basis,
                        n_gaussians=n_gaussians,
                        n_heads=n_heads,
                        cutoff_network=cutoff_network,
                        cutoff=cutoff,
                        activation=activation,
                        apply_transition_function=apply_transition_function,   
                        post_ln=post_ln,
                        resnet=resnet,
                    )
                    for _ in range(n_scales)    
                ]                       
        )
               
                       
        ### Time Embedding:
        self.time_embedding = _gen_timing_signal(n_interactions, n_atom_basis) ###  [1,n_interactions,n_atom_basis]
        
        self.time_embedding_list = torch.split(self.time_embedding,1,dim=1) ### n_interactions*[1,1,n_atom_basis]
            
        ### ACT:
        self.use_act = use_act            
        if self.use_act and n_interactions>=1 :                       
            self.act_blocks = nn.ModuleList(
                    [
                        AdaptiveComputationTime(
                            n_atom_basis=n_atom_basis,
                            activation=activation,
                            bias=1., ### -1., ### 1., ### Initial bias, positive so that no long-ponder in early stage.
                        )  
                        for _ in range(n_scales)    
                    ]                       
            )            
        
            
        #################
        # set attributes
        self.n_scales = n_scales
        self.use_act = use_act
        self.n_interactions = n_interactions           
        
        self.return_intermediate = return_intermediate
        self.charged_systems = charged_systems
        if charged_systems:
            self.charge = nn.Parameter(torch.Tensor(1, n_atom_basis))
            self.charge.data.normal_(0, 1.0 / n_atom_basis ** 0.5)

        self.out_net = nn.Sequential(
                MLP(n_atom_basis, 1, None, 2, swish),
            )
        self.pooling = SumPooling()

    def forward(self, g, cell=None):
        """Compute atomic representations/embeddings.

        Args:
            inputs (dict of torch.Tensor): SchNetPack dictionary of input tensors.

        Returns:
            torch.Tensor: atom-wise representation.
            list of torch.Tensor: intermediate atom-wise representations, if
            return_intermediate=True was used.

        """
        logger = logging.getLogger('rl.'+__name__)

        g.ndata['xyz'].requires_grad_()
        # get atom embeddings for the input atomic numbers
        e = self.embedding(g.ndata['z'].long())

        x = e

        x0 = x
        
        # store intermediate representations    
        xs = [x0]
                                   
        x_act_list = []  ### or [x]
        t_act_list = []
        
        # compute interaction block to update atomic embeddings       
        
        for i_scale in range(self.n_scales):        
            
            interaction_function = self.interaction_blocks[i_scale] 
            
            if self.use_act : 
                act_function = self.act_blocks[i_scale] 
            
            if self.use_act : 
                max_hop = self.n_interactions ### -1
                x, (remainders,n_updates) = act_function(g, x, e, interaction_function, self.time_embedding, max_hop, cell=cell)
           
            else:                    
                for i_interaction in range(self.n_interactions):                           
                    t = self.time_embedding_list[i_interaction]                
                    x = interaction_function(g, x, t, cell=cell) 
            
            ### After each Universal Transformer block:                                        
            xs.append(x)
            
        return x
                

