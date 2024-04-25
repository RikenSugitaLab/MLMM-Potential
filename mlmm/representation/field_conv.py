import torch
from torch import nn
from mlmm.utils.logger import setup_logger
from mlmm.representation.schnetconv import SchNetConv
from mlmm.nn.field_interaction.field_nn import  TensorInteraction, multipole_interaction
from mlmm.nn.acsf import  GaussianSmearing, LogGaussianDistribution 
from mlmm.nn.cfconv import distance
from mlmm.nn.graph_norm import gen_read_out, graph_norm
from mlmm.nn.cutoff import *
from mlmm.utils.utils import check_value, compute_interacton_tensor2
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
                 if_Tensor_Interaction=True,
                 activation="swish",
                 aggregator_mode:str='sum',
                 LODE = False,
                 accumulate_update=False,
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
        self.interactions_v = nn.ModuleList()
        self.update_v = nn.ModuleList()
        self.if_Tensor_Interaction = if_Tensor_Interaction
        self.dipole_module_list = nn.ModuleList()
        self.LODE = LODE
        self.gamma = gamma
        self.delta = delta
        self.accumulate_update = True

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
        
        for _ in range(self.n_interactions):
            self.interactions_v.append(message_v(n_basis,n_atom_basis,activation=activation,cutoff_network=cutoff_network,cross_product=cross_product))
            self.update_v.append(gated_update_v(n_atom_basis,activation=activation))
            self.dipole_module_list.append(gated_update_v(n_atom_basis,activation=activation))


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
            self.virtual_charge = nn.ModuleList([
                MLP(n_atom_basis, 1, activation=activation) for _ in range(n_interactions)
            ])

    def forward(self, g_orig, cell=None):
        logger = logging.getLogger('rl.'+__name__)
        try:
            g = dgl.edge_type_subgraph(g_orig.local_var(),[('qm','Hqm','qm')])
            g.set_batch_num_edges(g_orig.batch_num_edges(etype='Hqm'))
            g.set_batch_num_nodes(g_orig.batch_num_nodes(ntype='qm'))

            g2 = dgl.edge_type_subgraph(g_orig.local_var(),[('qm','Hqm2','qm')])
            g2.set_batch_num_edges(g_orig.batch_num_edges(etype='Hqm2'))
            g2.set_batch_num_nodes(g_orig.batch_num_nodes(ntype='qm'))
        except:
            g = g_orig.local_var()
            g2 = g_orig.local_var()

        ## calculate distance and distance expansion
        rij_all, vij_all = self.distance(g2,cell=cell)
        rij, vij = self.distance(g,cell=cell)
        fij = self.distance_expansion(rij)

        feat = self.embedding(g.ndata['z'].long())#.double()

        if check_value(feat):
            logger.critical('value error: nan or inf')
            assert check_value(feat)==0 

        feat_v = torch.zeros([feat.shape[0],feat.shape[1],3],device=feat.device)
        
        if self.LODE:
            T = compute_interacton_tensor2(g2, rij_all, vij_all, gamma=self.gamma, delta = self.delta)
            g = update_T(g,T)

        for l in range(self.n_interactions):
            v_total = self.interactions[l](g, feat, fij, rij) 

            m_v = self.interactions_v[l](g, feat, fij, vij/(1+rij[:,None]), rij, feat_v = feat_v)
            v_total_v = m_v

            feat = feat + v_total
            feat_v = feat_v + v_total_v 


            q, mu = self.dipole_module_list[l](feat, feat_v)

            m_mul_s, m_mul_v = self.multipole_layers[l](g, q, mu)
            feat = feat + m_mul_s
            feat_v = feat_v + m_mul_v
            v_total_v = v_total_v + m_mul_v
            v_total = v_total + m_mul_s

            # Dipole interaction
            if self.if_Tensor_Interaction:
                v_dipoles, v_dipoles_v = self.dipole_interactions[l](
                    g, mu, rij, vij/(1+rij[:,None]), f_ij=fij, feat=feat 
                )

                feat = feat + v_dipoles
                feat_v = feat_v + v_dipoles_v
                v_total_v = v_total_v + v_dipoles_v
                v_total = v_total + v_dipoles

            else:
                v_dipoles = 0

            mg_s, mg_v = self.update_v[l](feat, feat_v)

            feat = feat + mg_s
            feat_v = feat_v + mg_v
            v_total_v = v_total_v + mg_v
            v_total = v_total + mg_s

            if self.do_graph_norm:
                feat = self.graph_norm[l](g,feat)

            if self.LODE:
                if not self.accumulate_update:
                    g = update_T(g,T,aug=False)
                T = compute_interacton_tensor2(g2, rij_all, vij_all, self.virtual_charge[l](v_total).squeeze(), gamma=self.gamma, delta=self.delta)
                g = update_T(g,T)
            
        if self.graph_rep:
            graph_rep = self.rep_pooling(g, feat)
        else:
            graph_rep = None

        return feat, feat_v, graph_rep
    
def update_T(g, T, aug=True):
    for key in key_list:
        if key in T.keys():
            if aug:
                g.ndata[key] = g.ndata[key] + T[key]
            else:
                g.ndata[key] = g.ndata[key] - T[key]
    return g
            


