# Manual for ML/MM Model 

## Model Part
wf: weight for loss term of force  
wq: weight for loss term of charge of qm atoms  
total_charge: the net charge of QM region. Default: 0.0
train_total_energy: whether the total energy/force of QM region will be used for training. If it's set to false, the difference of total energy/force between couloumb energy/force calculated from QM charge output from Quantum Calculation will be used for training. Default: True  
reference_model: Wherther couloumb energy/force calculated from predicted QM charge will be used in Readout Module. when this argument is set to true, then please store the energy and force predicted by reference model (Coulomb model using QM charge derived from QM calculation) into label file (.h5), the corresponding key is ref_energy and ref_force. Default: True.

### rep_model:
n_interactions: number of interaction layers. Default: 6  
n_atom_basis: dimension of invariant/equivariant feature vector. Default: 256  
n_filters: dimension of filtering space (commonly the same as n_atoms_basis. Default: 256  
n_basis: number of radial basis functions. Default: 128  
max_z: max nuclear charge of elements in the system. Default: 100  
cutoff: cutoff distance for graph convolution. Default: 3.0 angstrom  
basis_func: the radial basis function used, "gaussian" (Gaussian basis) and "sin" (Painn radial basis function). Default: gaussian   
cutoff_network: the type of cut off function. "MollifierCutoff" and "CosineCutoff" are available. Default: MollifierCutoff    
do_graph_norm: whether a graph normalization layer is used at the end each interction layer. Default: False   
if_Tensor_Interaction: If dipole-dipole message is used. Default: True  
LODE: whether long-range descriptor contributed from QM atoms are considered. Default: True  
gamma: parameter of switch function in LODE (Unit: Bohr)  
delta: parameter of switch function in LODE (Unit: Bohr)  
accumulate_update: Whether the update from QM atoms to long-range descriptor is accumulated. Default: False  
activation function: the activation function used in representation module(only support "swish"). Default: swish  

### out_model:
activation function: the activation function used in read out module (only support "shifted-plus" and "swish"). Default: swish  
n_layers: number of fully connected layers in MLP. Default: 2  

## Data Part
batch_size: batch size  
num_workers: the same option in dataloader (refer to pytorch)  
pin_memory: the same option in dataloader (refer to pytorch)  
persistent_worker: the same option in dataloader (refer to pytorch)  

# Dataset
datafile: path of trjectory file. Avaiable format: \*.filelist (a list of path of multiple trajectory files), all the format supported by MDAnalysis  
top: path of topology file. Format: All the format supported by MDAnalysis  
label_file: the path of label file. (hdf5)  
select_group: QM region. (Please refer to MDAnalysis about selection language)  
full_connected: whether all the qm atoms is connected with each other in graph. Default: False  
in_memory: whether all the trajectory data are loaded in memory. Default: False  
T_max_order: the truncation order of Talyor Expansion of external potential. Default: 3   
sr_raidus: only force of mm atoms within this distance (QM-MM distance) will be used for training. Default: 10 angstrom  

sharing_data: whether data server shared across dataloader processes are contructed. Default: True  
When num_workers is larger than 1 and multiple gpu cards are used, the data will replicated in each process, resulting in huge memory cost. To solve this problem, we utilize the shared_tensor function inf deep graph library to share the data accross the processes to reduce the memory cost. Moreover, to avoid the repeated calculation of Talyor Expansion of external potential and graph construction, we contruct a server to store these calculated data. After one epoch, those data need not to be calculated again.   

n_atom_mm: the max number of mm atoms within sr_raidius in all the data  
n_edge_qm_max: the max number of edges in QM region in all the data  
The two options are used to contruct data server.   

## ** Please only modify the options listed in this text and leave out other ones shown in config.yml, other options are either experimental ones or deprecated ** 
