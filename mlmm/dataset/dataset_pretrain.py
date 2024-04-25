import torch 
import logging
import dgl
import os
from dgl import save_graphs, load_graphs
import netCDF4 as nclib
from dgl.data.utils import makedirs, save_info, load_info
from dgl.data import DGLDataset
import MDAnalysis as mda
import numpy as np
import MDAnalysis.lib.nsgrid as ngrid
import time
from torch.utils.data import DataLoader
from MDAnalysis.coordinates.memory import MemoryReader
from mlmm.utils.utils import read_filelist
if "use_fp16" in os.environ:
    if os.environ['use_fp16'] == 'True':
        use_fp16=True
    elif os.environ['use_fp16'] == 'False':
        use_fp16=False
    else:
        raise AssertionError("wrong setting, use_fp16 can only be True or False")
else:
    use_fp16=False

__all__=["build_graph","my_collate"]

def build_graph(graph_data,
                topfile,
                datafile,
                idx_file=None,
                select_group='all',
                in_memory=True
                ):
    
    if not isinstance(graph_data,tuple) and graph_data is not None:
        raise Exception("graph_data should be a tuple",graph_data)

    format = datafile.split('.')[-1]
    datafile = datafile

    if format == 'filelist':
        datafile = read_filelist(datafile)
        if datafile[0].split('.')[-1]=='dcd':
            nc = mda.Universe(topfile,datafile,in_memory=in_memory)
        if datafile[0].split('.')[-1]=='nc':
            nc = mda.Universe(topfile,datafile,in_memory=in_memory)
            #f = nclib.MFDataset(datafile)
            #nc = mda.Universe(topfile,f.variables['coordinates'][:],format=MemoryReader)
            #dimensions_array = np.concatenate([f.variables['cell_lengths'][:],f.variables['cell_angles'][:]],axis=1)
            #nc.trajectory.dimensions_array = dimensions_array
            #f.close()
    else:
        nc = mda.Universe(topfile,datafile,in_memory=in_memory)


    if idx_file is not None:
        idxb = np.loadtxt(idx_file,dtype=np.int32)

    if graph_data is None:
        if idx_file is None:
            n_frame = nc.trajectory.n_frames
        else:
            n_frame = idxb.shape[0]
        graph_data = tuple((np.arange(n_frame),np.arange(n_frame)))

    select_group = nc.select_atoms(select_group)
    element={'F':9,'O':8,'C':6,'N':7,'S':16,'H':1,'Cl':35}
    atom_type= select_group.atoms.elements
    atom_index=np.array([element[i] for i in atom_type],dtype=np.int32)

    idx = select_group.ix
    g = dgl.graph(graph_data)

    if in_memory:
        atom_index = np.tile(atom_index[None,:],(nc.trajectory.n_frames,1))
        if idx_file is None:
            idxb = np.arange(nc.trajectory.n_frames)

        if use_fp16:
            g.ndata['xyz'] = torch.from_numpy(nc.trajectory.coordinate_array[idxb,idx]).half()
            g.ndata['cell'] = torch.from_numpy(nc.trajectory.dimensions_array[idxb]).half()
        else:
            g.ndata['xyz'] = torch.from_numpy(nc.trajectory.coordinate_array[idxb,idx])
            g.ndata['cell'] = torch.from_numpy(nc.trajectory.dimensions_array[idxb]).float()

        g.ndata['z'] = torch.from_numpy(atom_index[idxb])
    else:
        if idx_file is not None:
            trj_filter = nc.trajectory[idxb]
        else:
            trj_filter = nc.trajectory

        coor = []
        dimension = []
        atom_index2 = []
        for ts in trj_filter: 
            coor.append(nc.atoms.positions[idx])
            dimension.append(nc.atoms.dimensions)
            atom_index2.append(atom_index)

        coor = np.stack(coor)
        dimension = np.stack(dimension)
        atom_index2 = np.stack(atom_index2)

        if use_fp16:
            g.ndata['xyz'] = torch.from_numpy(coor).half()
            g.ndata['cell'] = torch.from_numpy(dimension).half()
        else:
            g.ndata['xyz'] = torch.from_numpy(coor)
            g.ndata['cell'] = torch.from_numpy(dimension).float()

        g.ndata['z'] = torch.from_numpy(atom_index2)
    return g, nc.trajectory.n_frames

def my_collate(input_nodes,output_nodes,blocks,cutoff=4):
    total_size1 = np.unique(np.concatenate([input_nodes,output_nodes])).shape[0]
    total_size2 = input_nodes.shape[0]

    assert total_size1 == total_size2, "error: input nodes don't contain output nodes"
    assert (input_nodes[:output_nodes.shape[0]] - output_nodes).sum() == 0, "error: output nodes don't lie on the top of input nodes"

    n_atoms = blocks[0].srcdata['xyz'].shape[1]

    graph_list = []
    cell_list = []
    for i in range(total_size2):
        neighbor_ob=ngrid.FastNS(cutoff,blocks[0].srcdata['xyz'][i].numpy(),blocks[0].srcdata['cell'][i].numpy())
        ns_result = neighbor_ob.self_search()
        neighbor_list=ns_result.get_pairs()
        dis = ns_result.get_pair_distances()
        neighbor_list=neighbor_ob.self_search().get_pairs()
        g=dgl.graph((torch.IntTensor(neighbor_list[:,0]),torch.IntTensor(neighbor_list[:,1])),num_nodes=n_atoms)
        g.ndata['xyz'] = blocks[0].srcdata['xyz'][i] 
        g.ndata['z'] = blocks[0].srcdata['z'][i] 
        g.edata['dis']=torch.from_numpy(dis)
        graph_list.append(g)

    return graph_list, blocks[0].srcdata['cell']
    


        
        

