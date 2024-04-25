import torch
from torch import nn
from mlmm.utils.logger import setup_logger
from mlmm.representation.schnetconv import SchNetConv
from mlmm.nn.field_interaction.field_nn import  Dipole_Layer, TensorInteraction, multipole_interaction
from mlmm.nn.acsf import GTO_Expansion, GaussianSmearing, LogGaussianDistribution 
from mlmm.nn.cfconv import distance
from mlmm.nn.graph_norm import gen_read_out, graph_norm
from mlmm.nn.cutoff import *
from mlmm.utils.utils import check_value, compute_interacton_tensor, compute_interacton_tensor2
from mlmm.nn.aggregator import aggregator
import logging
import dgl
from mlmm.nn.blocks import MLP
from mlmm.nn.activations import swish,swish2
from mlmm.nn.update_coor import Update_coor
from mlmm.nn.equvariant_message import *

__all__ = ["FieldNet"]
key_list = ['T0','T1','T2','T3']

class FieldNet(nn.Module):
    """SchNet architecture for learning representations of atomistic systems.

    Args:
        n_interactions (int, optional): number of interaction blocks.
        n_atom_basis (int, optional): number of features to describe atomic environments.
            This determines the size of each embedding vector; i.e. embeddings_dim.
        n_filters (int, optional): number of features in continuous filter
        n_basis (int, optional): number of basis functions used to expand
            atomic distances.
        dipole_features (int, optional): number of features used to represent dipole).
        max_z (int, optional): maximum nuclear charge allowed in database. This
            determines the size of the dictionary of embedding; i.e. num_embeddings.
        cutoff (float, optional): cutoff radius.
        L (int): parameter for GTO expansion
        basis_func (str, optional): the basis function used to expand distance
        cutoff_network (str, optional): cutoff layer.
        dipole_cutoff (nn.Module, optional): dipole cutoff layer.
        do_graph_norm (bool, optional): if do graph normalization.
        graph_rep (bool, optional): if get graph representation.
        elec_charge(bool, optional): if adding charge-field message.
    """

    def __init__(self,
                 n_interactions=3,
                 n_atom_basis=100,
                 n_filters=100,
                 n_basis=25,
                 dipole_features=100,
                 max_z=100, 
                 cutoff=5,
                 basis_func: str = "gaussian",
                 cutoff_network='MollifierCutoff',
                 dipole_cutoff=None,
                 do_graph_norm=False,
                 graph_rep=False,
                 cross_product = False,
                 equvariant_message = False,
                 if_Tensor_Interaction=True,
                 activation="swish",
                 aggregator_mode:str='sum',
                 LODE = False,
                 gamma = 0.37,
                 delta = 0.5
                 ):
        super().__init__()
        if activation=="swish":
            activation = swish
        elif activation=="swish2":
            activation = swish2

        self.n_interactions = n_interactions

        self.n_atom_basis = n_atom_basis
        self.cutoff = cutoff
        self.embedding = nn.Embedding(max_z, n_atom_basis, padding_idx=0)
        self.interactions = nn.ModuleList()
        self.do_graph_norm = do_graph_norm
        self.cross_product = cross_product
        self.field_evo = nn.ModuleList()
        self.interactions_v = nn.ModuleList()
        self.update_v = nn.ModuleList()
        self.if_Tensor_Interaction = if_Tensor_Interaction
        self.equvariant_message = equvariant_message
        self.dipole_module_list = nn.ModuleList()
        self.LODE = LODE
        self.gamma = gamma
        self.delta = delta

        if do_graph_norm:
            self.graph_norm = nn.ModuleList() 

        self.distance = distance()

        if dipole_cutoff is None:
            dipole_cutoff = cutoff

        if cutoff_network == "MollifierCutoff":
            cutoff_network = MollifierCutoff(cutoff = cutoff)
            dipole_cutoff_network = MollifierCutoff(cutoff = dipole_cutoff)
        elif cutoff_network == "CosineCutoff":
            cutoff_network = CosineCutoff(cutoff = cutoff)
            dipole_cutoff_network = CosineCutoff(cutoff = dipole_cutoff)

        self.basis_func = basis_func
        if basis_func=="loggausssian":
            self.distance_expansion = LogGaussianDistribution(d_max=cutoff,num_rbf=n_basis)
        elif basis_func == "gaussian":
            self.distance_expansion = GaussianSmearing(
                0.0, cutoff, n_gaussians=n_basis)
        
        if equvariant_message:
            for _ in range(self.n_interactions):
                self.interactions_v.append(message_v(n_basis,n_atom_basis,activation=activation,cutoff_network=cutoff_network))
                self.update_v.append(gated_update_v(n_atom_basis,activation=activation))
                self.dipole_module_list.append(dipole_layer_eq(n_atom_basis,activation=activation))


        for _ in range(self.n_interactions):
            self.interactions.append(SchNetConv(
                                          n_atom_basis=n_atom_basis,
                                          n_filters=n_filters,
                                          n_basis=n_basis,
                                          cutoff_network=cutoff_network,
                                          activation=activation,
                                          aggregator=aggregator(mode=aggregator_mode)
                                          ))
            if do_graph_norm:
                self.graph_norm.append(graph_norm(n_atom_basis))

        if not equvariant_message:
            self.dipole_update = nn.ModuleList([
                                 Dipole_Layer(n_atom_basis, dipole_features,
                                 transform=True,
                                 activation=activation,
                                 cutoff_network=dipole_cutoff_network, cross_product=cross_product)
                                 for _ in range(n_interactions + 1)
                                    ])

        if if_Tensor_Interaction:
            self.dipole_interactions = nn.ModuleList([
                    TensorInteraction(dipole_features,n_atom_basis,cutoff=cutoff, activation=activation,
                                        n_basis=n_basis, aggregator = aggregator(mode=aggregator_mode)) for _ in range(n_interactions)
            ])

        self.multipole_layers = nn.ModuleList([
            multipole_interaction(dipole_features, n_atom_basis, activation=activation) for _ in range(n_interactions)
        ])
        self.graph_rep = graph_rep
        if graph_rep:
            self.rep_pooling = gen_read_out(agg='softmax')
        
        if self.LODE:
            self.virtual_charge = MLP(n_atom_basis,1,activation=activation)

    def forward(self, g_orig, cell=None):
        logger = logging.getLogger('rl.'+__name__)
        g = g_orig.local_var()

        ## calculate distance and distance expansion
        rij, vij = self.distance(g,cell=cell)

        fij = self.distance_expansion(rij)

        # if (rij>self.cutoff+0.1).sum()>0:
            # logger.debug(f"{(rij>self.cutoff).sum()}")
            # logger.debug(f"{rij[rij>self.cutoff]}")
            # logger.critical('error: distance calculation errror, distance should not be large than cutoff distance')
            #assert (rij>self.cutoff).sum()==0

        feat = self.embedding(g.ndata['z'].long())#.double()

        if check_value(feat):
            logger.critical('value error: nan or inf')
            assert check_value(feat)==0 

        if not self.equvariant_message:
            mu0 = self.dipole_update[0](g, feat, rij, vij/(1+rij[:,None]))
            mu = mu0
        #mu = self.graph_norm_dipole[0](g,mu)
            if check_value(mu):
                logger.critical('value error: nan or inf')
                assert check_value(mu)==0 
        else:
            # feat_v = self.interactions_v[0]()
            feat_v = torch.zeros([feat.shape[0],feat.shape[1],3],device=feat.device)
        
        if self.LODE:
            T = compute_interacton_tensor2(g, rij, vij, gamma=self.gamma, delta = self.delta)
            g = update_T(g,T)

        for l in range(self.n_interactions):
            v_total = self.interactions[l](g, feat, fij, rij) 
            print(f"{l} interaction-v_total:{abs(v_total).max()}")

            if check_value(v_total):
                logger.critical('value error: nan or inf')
                assert check_value(v_total)==0 
            
            if self.equvariant_message:
                m_v = self.interactions_v[l](g, feat, fij, vij/(1+rij[:,None]), rij, feat_v = feat_v)
                feat_v = feat_v + m_v
                v_total_v = m_v
                print(f"{l} interaction-m_v:{abs(m_v).max()}")

                if check_value(m_v):
                    logger.critical('value error: nan or inf')
                    assert check_value(m_v)==0 

                mu = self.dipole_module_list[l](feat, feat_v)
                print(f"{l} interaction-mu:{abs(mu).max()}")
                if check_value(mu):
                    logger.critical('value error: nan or inf')
                    assert check_value(mu)==0 

            m_mul_s, m_mul_v = self.multipole_layers[l](g, feat, mu)
            v_total_v = m_mul_v
            print(f"{l} interaction-m_mul_s:{abs(m_mul_s).max()}")
            print(f"{l} interaction-m_mul_v:{abs(m_mul_v).max()}")

            if check_value(m_mul_s) or check_value(m_mul_v):
                logger.critical('value error: nan or inf')
                assert check_value(m_mul_s)==0 
                assert check_value(m_mul_v)==0 
        
            # Dipole interaction
            if self.if_Tensor_Interaction:
                v_dipoles, v_dipoles_v = self.dipole_interactions[l](
                    g, mu, rij, vij/(1+rij[:,None]), f_ij=fij, feat=feat 
                )
                # v_dipoles, v_dipoles_v = self.dipole_interactions[l](
                    # g, mu, rij, vij, f_ij=fij, feat=feat 
                # )
                print(f"{l} interaction-v_dipoles:{abs(v_dipoles).max()}")
                print(f"{l} interaction-v_dipoles_v:{abs(v_dipoles_v).max()}")

                v_total_v  = v_total_v + v_dipoles_v
                if check_value(v_dipoles): 
                    logger.critical('value error: nan or inf')
                    assert check_value(v_dipoles)==0 
            else:
                v_dipoles = 0

            # Update features
            v_total = v_total + m_mul_s + v_dipoles
            if self.equvariant_message:
                feat = feat + v_total
                feat_v = feat_v + v_total_v 

                mg_s, mg_v = self.update_v[l](feat, feat_v)
                print(f"{l} interaction-mg_v:{abs(mg_v).max()}")
                print(f"{l} interaction-mg_s:{abs(mg_s).max()}")

                if check_value(mg_s) or check_value(mg_v):
                    logger.critical('value error: nan or inf')
                    assert check_value(mg_s)==0 
                    assert check_value(mg_v)==0 

                feat = feat + mg_s
                feat_v = feat_v + mg_v
            else:
                feat = feat + v_total

            if self.do_graph_norm:
                feat = self.graph_norm[l](g,feat)

            if not self.equvariant_message:
                mu_update = self.dipole_update[l + 1](g, v_total, rij, vij/(1+rij[:,None]))
                mu = mu + mu_update

            if self.LODE:
                g = update_T(g,T,aug=False)
                T = compute_interacton_tensor2(g, rij, vij, self.virtual_charge(v_total).squeeze(), gamma=self.gamma, delta=self.delta)
                g = update_T(g,T)

            if check_value(mu):
                logger.critical('value error: nan or inf')
                assert check_value(mu)==0 
            
        if self.graph_rep:
            graph_rep = self.rep_pooling(g, feat)
        else:
            graph_rep = None

        return feat, feat_v, graph_rep, mu, rij
    
def update_T(g, T, aug=True):
    for key in key_list:
        if key in T.keys():
            if aug:
                g.ndata[key] = g.ndata[key] + T[key]
            else:
                g.ndata[key] = g.ndata[key] - T[key]
    return g
            


