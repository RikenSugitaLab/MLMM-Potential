from collections import defaultdict
from multiprocessing.sharedctypes import RawArray
import pickle
from multiprocessing import synchronize
from multiprocessing.connection import Listener
import torch as th
import os
import dgl
from dgl.multiprocessing.pytorch import shared_tensor
import netCDF4 as nc
from dgl.data.utils import makedirs, save_info, load_info
from dgl.data import DGLDataset
import MDAnalysis as mda
import numpy as np
import MDAnalysis.lib.nsgrid as ngrid
import multiprocessing as mp
from torch.utils.data import Dataset, DataLoader
from MDAnalysis.coordinates.memory import MemoryReader
from mlmm.utils.utils import read_filelist, write_data
from copy import deepcopy
import glob
import h5py
import re
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
from memory_profiler import profile
import psutil

def get_current_memory_gb():
# 获取当前进程内存占用。
    pid = os.getpid()
    p = psutil.Process(pid)
    info = p.memory_full_info()
    return info.uss / 1024. / 1024. / 1024.
# from mlmm.baal.get_label import DftbCalculator

qminfo2genesis = 0.280028557027
CONV_UNIT_ENE = 627.5095
CONV_UNIT_LEN    = 0.52917721092
CONV_UNIT_FORCE  = CONV_UNIT_ENE / CONV_UNIT_LEN


__all__=["collate","mydataset"]

def collate(samples):

    if len(samples[0])>2:
        graphs, keyname, label, cell = map(list,zip(*samples))
        keyname = keyname[0]
    else:
        graphs, cell = map(list,zip(*samples))
        label=None
            
    batched_graph = dgl.batch(graphs)
    cell = th.stack(cell)

    labelout = {}
    if len(samples[0])>2:
        label = list(zip(*label))

        for i in range(len(label)):
            labelout[keyname[i]] = th.stack(label[i])

    return batched_graph, labelout, cell 


### building the dataloader (read data from file, provide data to model)
class mydataset(Dataset):#DGLDataset):
    def __init__(self,
                 datafile = '',         # trajectory file   
                 top='',                # topology file
                 label_file='',         # label file (such as energy charge or external field)
                 select_group='all',    # the atom groups which is used as input
                 full_connected=False,  # whether all atoms are connected by an edge
                 cutoff=4,              # cutoff distance to set an edge
                 in_memory=True,        # whether transfer the data into memory cache
                 force_field_type=False,# whether embedding atoms by it's type in force field
                 add_self_loop=False,   # whether self loop edge is added
                 if_compute_elec = True, 
                 transform_unit = True, 
                 columb_energy = False,
                 sr_raidus = 10.0, 
                 sharing_data = True,
                 T_max_order = 1,
                 analytical_grad = True,
                 n_atom_mm_max = 0,
                 n_edge_qm_max = 0,
                 server_ckpt = None
                 ):
        """molecule dataset for predicting properties based on MDAnalysis 
           supported format for datafile: netcdf, dcd, npz

        Args:
            datafile (str) : path of trajectroy file
            top (str): path of topology file
            label_file (str): path of label file
            select_group (str): the atom group used as input
            full_connected (bool): whether all atoms are connected 
            cutoff (float): cutoff distance to set an edge
            in_memory (bool): whether transfer the data into memory cache
            force_field_type (bool): whether embedding atoms by it's type in force field
            add_self_loop (bool): whether self loop edge is added
            if_compute_elec (bool): whether compute external electric field and strength during traing

        """
        self.datafile = datafile
        self.cutoff = cutoff
        self.read_vel = False
        self.g = []
        self.force_field_type = force_field_type
        self.add_self_loop = add_self_loop
        self.if_compute_elec = if_compute_elec 
        self.select_index = select_group
        self.num_classes = 25 ### need modified
        self.transform_unit = transform_unit
        self.columb_energy = columb_energy
        self.sr_raidus = sr_raidus
        self.T_max_order = T_max_order
        self.analytical_grad = analytical_grad
        self.in_memory = in_memory
        self.top = top
        self.select_group = select_group
        self.full_connected = full_connected
        self.label_file = label_file
        self.sharing_data = sharing_data
        self.n_atoms_mm_max = n_atom_mm_max
        self.n_edge_qm_max = n_edge_qm_max
        self.server_ckpt = server_ckpt
        self.ckpt_path = None

        self.trj_format = ['dcd','nc','filelist','npz']

        ### reading the trajectory file ###
        if isinstance(datafile, str):
            self.format = datafile.split('.')[-1]
        else:
            self.format = ''

    # @profile
    def load_data(self):
        print('loading the data')
        self.ddp = th.distributed.is_initialized()
        
        if self.ddp:
            self.main_rank = (th.distributed.get_rank()==0)
        else:
            self.main_rank = True

        if self.format in self.trj_format and not isinstance(self.datafile, mda.Universe):

            if self.format == 'filelist':
               self.datafile = read_filelist(self.datafile)

            if self.format != 'npz':
                print('loading the nc')
                self.nc = mda.Universe(self.top,self.datafile,in_memory=self.in_memory,forces=True)
                if self.nc.atoms.ids[0] != 0:
                    self.nc.atoms.ids = self.nc.atoms.ids - 1
                print('done')
                self.n_frames = self.nc.trajectory.n_frames

                if hasattr(self.nc,'velocity_array'):
                    if self.nc.trajectory.velocity_array is not None:
                        self.read_vel=True
            elif self.format == 'npz':
                data = np.load(self.datafile)
                if self.top != '':
                    self.nc = mda.Universe(self.top,data['R'],format=MemoryReader)
                else:
                    n_atoms = data['R'].shape[1]
                    self.nc = mda.Universe.empty(n_atoms,trajectory=True,velocities=True,forces=True)
                    self.nc.load_new(data['R'],format=MemoryReader)
                    self.atom_index = data['z'][0]

                self.n_frames = data['R'].shape[0]
                if 'vel' in data:
                    self.nc.trajectory.velocity_array=data['vel']
                    self.read_vel = True

                if 'dimensions' in data:
                    self.nc.dimensions=data['dimensions']
        elif isinstance(self.datafile, mda.Universe):
            self.nc = self.datafile
            self.n_frames = self.nc.trajectory.n_frames
        else:
            raise NotImplementedError('data format not supported: {}',format(self.format))

        self.select_group = self.nc.select_atoms(self.select_group)
        self.qmidx = self.select_group.ix_array
        self.mmidx = np.setdiff1d(np.arange(self.nc.atoms.n_atoms), self.qmidx)
        self.n_atoms_mm = self.mmidx.shape[0]
        self.n_atoms_qm = self.qmidx.shape[0]
        if self.top != '' or isinstance(self.datafile, mda.Universe):
            if self.force_field_type: 
                key_name = np.unique(self.select_group.atoms.names)
                element={}
                for i in range(key_name.shape[0]):
                   element[key_name[i]]=i
                atom_type= self.nc.atoms.names
            else:
                element={'F':9,'O':8,'C':6,'N':7,'S':16,'H':1,'CL':35, 'LA':88,'P':15}
                element2=['CL','LA']
                atom_type= self.select_group.atoms.names

            atom_index = []
            for i in atom_type:
                if i[:2] in element2:
                    atom_index.append(element[i[:2]])
                else:
                    atom_index.append(element[i[0]])
            self.atom_index = np.array(atom_index,dtype=np.float32)
            self.atom_mass=np.array(self.select_group.atoms.masses)

        print('loading the h5py')
        if self.label_file != '' and isinstance(self.label_file,str):
            if self.label_file.split('.')[-1]=='npz':
               label_data = np.load(self.label_file)
               self.label_data = {}
               for key in list(label_data):
                   if key in ['z','R','v']:
                       continue
                   self.label_data[key] = th.from_numpy(label_data[key])
                   assert self.label_data[key].shape[0]==self.n_frames, 'shape mismatch'
            elif self.label_file.split('.')[-1]=='txt':
               self.label_data =th.from_numpy(np.loadtxt(self.label_file))
               assert self.label_data.shape[0]==self.n_frames, 'shape mismatch'
            elif self.label_file.split('.')[-1]=='h5py' or self.label_file.split('.')[-1]== 'h5':
                f = h5py.File(self.label_file,'r')
                self.label_data = {}
                for key in f.keys(): 
                   if key in ['z','R','v']:
                       continue
                #    self.label_data[key] = th.from_numpy(f.get(key)[:][1:]).float()
                   if key == "extfield":
                       if self.sharing_data and self.ddp:
                            if self.ddp:
                                th.distributed.barrier()
                            self.label_data['T1'] = shared_tensor((self.n_frames,self.n_atoms_qm,3))
                            if self.main_rank:
                                #self.label_data['T1'] = th.from_numpy(f.get(key)[:]).float()
                                self.label_data['T1'].copy_(th.from_numpy(f.get(key)[:]).float())
                       else:
                           self.label_data['T1'] = th.from_numpy(f.get(key)[:]).float()
                       key = 'T1'
                   elif key == "extfield_st":
                       if self.sharing_data and self.ddp:
                            if self.ddp:
                                th.distributed.barrier()
                            self.label_data['T0'] = shared_tensor((self.n_frames,self.n_atoms_qm))
                            if self.main_rank:
                                self.label_data['T0'].copy_(th.from_numpy(f.get(key)[:]).float())
                       else:
                           self.label_data['T0'] = th.from_numpy(f.get(key)[:]).float()
                       key = 'T0'
                   else:
                       if self.nc.coord.has_forces and (key == "mm_force" or key == "qm_force"):
                          continue

                       if self.sharing_data and self.ddp and key != "mm_force":
                            if self.ddp:
                                th.distributed.barrier()
                            self.label_data[key] = shared_tensor(f[key].shape)
                            if self.main_rank:
                                self.label_data[key].copy_(th.from_numpy(f.get(key)[:]).float())
                       else:
                           self.label_data[key] = th.from_numpy(f.get(key)[:]).float()
                   assert self.label_data[key].shape[0]==self.n_frames, 'shape mismatch'
            else:
               raise Exception('can not process the label file ended with'+self.label_file.split('.')[-1])
        elif not isinstance(self.label_file, str):
            self.label_data = self.label_file
        else:
            self.label_data = None
        print('done')

        if isinstance(self.label_data, dict) and self.transform_unit:
            if self.main_rank or (not self.sharing_data):
                if 'mm_force' in self.label_data.keys():
                    self.label_data['mm_force'].copy_(self.label_data['mm_force']/CONV_UNIT_FORCE)

                if 'qm_force' in self.label_data.keys():
                    self.label_data['qm_force'].copy_(self.label_data['qm_force']/CONV_UNIT_FORCE)

                if 'qm_energy' in self.label_data.keys():
                    self.label_data['qm_energy'].copy_(self.label_data['qm_energy']/CONV_UNIT_ENE)

                if 'ref_energy' in self.label_data.keys():
                    self.label_data['ref_energy'].copy_(self.label_data['ref_energy']/CONV_UNIT_ENE)
                else:
                    self.label_data['ref_energy']=(th.zeros_like(self.label_data['qm_energy']).float())
                    self.label_data['ref_energy'].copy_(self.label_data['ref_energy']/CONV_UNIT_ENE)
                    Warning('no reference energy, please notice qm_energy should be compatible with option train_total_energy')

                if 'ref_force' in self.label_data.keys():
                    self.label_data['ref_force'].copy_(self.label_data['ref_force']/CONV_UNIT_FORCE)
                else:
                    self.label_data['ref_force']=(th.zeros([self.n_frames, self.n_atoms_qm,3]).float())
                    self.label_data['ref_force'].copy_(self.label_data['ref_force']/CONV_UNIT_FORCE)
                    Warning('no reference force, please notice qm_force should be compatible with option train_total_energy')

        if isinstance(self.label_data, dict) and self.transform_unit:
            if self.main_rank or (not self.sharing_data):
                if 'ref_energy' in self.label_data.keys():
                    self.label_data['qm_energy'].copy_(self.label_data['qm_energy'] - self.label_data['ref_energy'])

        self.get_statistic() 

        self.qmidx = np.sort(self.qmidx)
        self.mmidx = np.sort(self.mmidx)
        self.internal_idx = np.zeros(self.nc.atoms.n_atoms,dtype=np.int64)
        self.internal_idx[self.qmidx] = np.arange(self.qmidx.shape[0],dtype=np.int64) 
        self.internal_idx[self.mmidx] = np.arange(self.mmidx.shape[0],dtype=np.int64) 

        print('set up graph server')
        if self.sharing_data:
            g,cell,eid, diffidx = self.nglist_search(0,
                                            self.select_group,
                                            atom_index=self.atom_index,
                                            atom_mass=self.atom_mass,
                                            cutoff=self.cutoff,
                                            full_conneted=self.full_connected)
            T = self.compute_elec()
            for key in T.keys():
                if "Tij" not in key:
                    if len(g.ntypes)<1:
                        g.ndata[key] = T[key]
                    else:
                        g.ndata[key] = {'qm': th.from_numpy(T[key]).float()}
                else:
                    if self.analytical_grad:
                        g.edata[key] = {'Hqmmm':th.from_numpy(T[key][eid[0],eid[1]]).float()}
            self.graph_server = graph_server(self.main_rank,self.n_frames,self.n_atoms_mm_max,self.n_edge_qm_max,self.n_atoms_qm)
            self.graph_server.initialize(g)
            if self.server_ckpt is not None:
                server_state_dict = th.load(self.server_ckpt)
                self.graph_server.update(server_state_dict)
        print('done')
        
        super(mydataset,self).__init__()
                
    def get_statistic(self):
        if not self.columb_energy:
            self.mean_e = self.label_data['qm_energy'].mean()
            self.std_e = self.label_data['qm_energy'].std()
        else: 
            coor_qm = self.nc.trajectory.coordinate_array[:,self.qmidx]
            rij = np.linalg.norm(coor_qm[:,None] - coor_qm[:,:,None], axis=-1)
            cij = self.label_data['qm_charge'][:,None]*self.label_data['qm_charge'][:,:,None]
            columb_energy = (CONV_UNIT_LEN*cij.squeeze()/rij)
            columb_energy[np.isinf(columb_energy)] =  0.0
            columb_energy = (columb_energy.sum(-1).sum(-1))/2
            energy = self.label_data['qm_energy'] -  columb_energy
            self.mean_e = energy.mean()
            self.std_e = energy.std()

    # @profile
    def __getitem__(self,idx):
        self.nc.trajectory[idx]

        ### the frames are transformed into graph data on the fly when training ###

        if self.ddp and self.sharing_data:
            g = self.graph_server.get_graph(idx)
        else:
            g = None

        if g is None:
            g,cell,eid, diffidx = self.nglist_search(idx,
                                  self.select_group,
                                  atom_index=self.atom_index,
                                  atom_mass=self.atom_mass,
                                  cutoff=self.cutoff,
                                  full_conneted=self.full_connected)

            if self.add_self_loop:
                g = g.add_self_loop()

            if self.analytical_grad:
                T = self.compute_elec()
            else:
                T = self.compute_elec(mmidx=diffidx)

            for key in T.keys():
                if "Tij" not in key:
                    if len(g.ntypes)<1:
                        g.ndata[key] = T[key]
                    else:
                        g.ndata[key] = {'qm': th.from_numpy(T[key]).float()}
                else:
                    if self.analytical_grad:
                        g.edata[key] = {'Hqmmm':th.from_numpy(T[key][eid[0],eid[1]]).float()}
            
            if self.ddp and self.sharing_data:
                self.graph_server.append(g,idx)
            # print(f"idx:{idx}, rank:{th.distributed.get_rank()},have not g \n")
            # print(f"idx:{idx}, rank:{th.distributed.get_rank()}, mp_rank:{curPorc._identity[0]} have not g \n")
        else: 
            # print(f"idx:{idx}, rank:{th.distributed.get_rank()},get g \n")
            if self.nc.atoms.dimensions is None:
                cell = th.Tensor([999.0,999.0,999.0,90.0,90.0,90.0])
            else:
                cell = self.nc.atoms.dimensions

        if self.label_data is not None:
            label = {}
            keyname = []
            for key in self.label_data:
                if key not in ['R','v','T0','T1','T2','T3','T4']:

                   if key == 'qm_force':
                      if self.label_data[key][idx].ndim != 2:
                          a = 3
                   elif key == 'qm_energy':
                      if self.label_data[key][idx].ndim != 0:
                          a = 3
                   elif key == 'qm_charge':
                      if self.label_data[key][idx].ndim != 2:
                          a = 3
                   label[key] = self.label_data[key][idx]
                   keyname.append(key)
        
        if self.sharing_data:
            if self.graph_server.done and self.graph_server.save_time < 1:
                if self.ckpt_path is not None:
                    filename = os.path.join(self.ckpt_path,'graph_server.ckpt')
                else:
                    filename = './graph_server.ckpt'
                self.graph_server.save(filename)

        if self.label_data is not None:
            return g, keyname, tuple(label.values()), cell
        else:
            return g, cell 

    def __len__(self):
        return self.n_frames
    
    def nglist_search(self,frame_idx,atom_groups,atom_index,cutoff,atom_mass=None,full_conneted=False):
        ### transform the frames into graph data using mdanalysis ###
        n_atoms = atom_groups.n_atoms
        # atom_groups.universe.trajectory[frame_idx]

        if not full_conneted:
            if atom_groups.dimensions is None:
                cell = th.Tensor([999.0,999.0,999.0,90.0,90.0,90.0])
                Warning('no periodic box information, use default cell')
        ### neighbor list ###
            else:
                if atom_groups.dimensions[:3].sum()==0:
                    cell = th.Tensor([999.0,999.0,999.0,90.0,90.0,90.0])
                else:
                    cell = th.from_numpy(atom_groups.dimensions)

            neighbor_ob=ngrid.FastNS(cutoff,atom_groups.positions,cell)
            neighbor_list=neighbor_ob.self_search().get_pairs()
            snode = np.concatenate([neighbor_list[:,0],neighbor_list[:,1]])
            tnode = np.concatenate([neighbor_list[:,1],neighbor_list[:,0]])
        else:
            cell = th.Tensor([999.0,999.0,999.0,90.0,90.0,90.0])

            ### construct graph ###
            tnode = np.array([range(n_atoms)] * n_atoms)
            tnode = tnode[~np.eye(tnode.shape[0], dtype=bool)].reshape(
                                tnode.shape[0], -1
                               )
            tnode = th.from_numpy(tnode).reshape(-1)
            snode = (np.repeat(np.arange(n_atoms)[:,None],n_atoms-1,axis=1)).reshape(-1)


        if self.sr_raidus > 0.0:
            if not full_conneted:
                tnode2 = np.array([range(n_atoms)] * n_atoms)
                tnode2 = tnode2[~np.eye(tnode2.shape[0], dtype=bool)].reshape(
                                tnode2.shape[0], -1
                               )
                tnode2 = th.from_numpy(tnode2).reshape(-1)
                snode2 = (np.repeat(np.arange(n_atoms)[:,None],n_atoms-1,axis=1)).reshape(-1)
            else:
                tnode2 = tnode
                snode2 = snode
            sr_neighbor_ob = AtomNeighborSearch(self.nc.atoms[self.mmidx],box = cell)
            sr_atom = sr_neighbor_ob.search(atom_groups,self.sr_raidus)
            n_atoms_sr = sr_atom.n_atoms

            tnode_qm = (th.arange(n_atoms)[:,None]).repeat(1,n_atoms_sr).reshape(-1)
            snode_mm = ((th.arange(n_atoms_sr)[:,None]).repeat(1,n_atoms).T).reshape(-1)
            data_dict = {
                            ('qm', 'Hqm', 'qm'):  (snode, tnode),
                            ('qm', 'Hqm2', 'qm'): (snode2, tnode2),
                            ('mm', 'Hqmmm', 'qm'): (snode_mm,tnode_qm),
                        }
            ids = (self.internal_idx[self.qmidx[tnode_qm]],self.internal_idx[sr_atom.atoms.ids[snode_mm]])
            g=dgl.heterograph(data_dict)
            g.ndata['c'] = {'mm':th.from_numpy(np.float32(sr_atom.charges)),'qm':th.from_numpy(np.float32(atom_groups.charges))}
            g.ndata['xyz'] = {'mm': th.from_numpy(sr_atom.positions), 'qm': th.from_numpy(atom_groups.positions)}

            if hasattr(sr_atom,"forces"):
                if self.transform_unit:
                    g.ndata['f'] = {'mm': th.from_numpy(sr_atom.forces)/CONV_UNIT_FORCE, 'qm': th.from_numpy(atom_groups.forces)/CONV_UNIT_FORCE}
                else:
                    g.ndata['f'] = {'mm': th.from_numpy(sr_atom.forces), 'qm': th.from_numpy(atom_groups.forces)}
            else:
                mm_force = self.label_data['mm_force'][frame_idx,self.internal_idx[sr_atom.atoms.ids]]
                g.ndata['f'] = {'mm': mm_force, 'qm': self.label_data['qm_force'][frame_idx]}

            if hasattr(sr_atom,"velocities"):
                g.ndata['v'] = {'mm': th.from_numpy(sr_atom.velocities), 'qm': th.from_numpy(atom_groups.velocities)}

            g.ndata['z'] = {'qm': th.from_numpy(np.int32(atom_index))}
            g.ndata['m'] = {'qm': th.from_numpy(atom_mass)}
        else:
            g=dgl.graph((snode,tnode))
            ### node feature ###
            g.ndata['xyz'] = th.from_numpy(atom_groups.coordinates)

            if hasattr(sr_atom,"velocities"):
                g.ndata['v'] = th.from_numpy(atom_groups.velocities)
            g.ndata['z'] = th.from_numpy(np.int32(atom_index))

            if hasattr(sr_atom,"forces"):
                g.ndata['f'] = th.from_numpy(atom_groups.forces)/CONV_UNIT_FORCE
            else:
                g.ndata['f'] = self.label_data['qm_force']
            g.ndata['m'] = th.from_numpy(atom_mass)
            ids = None

            #cell = th.from_numpy(cell)

        # if charge is not None:
            # g.ndata['charge'] = th.from_numpy(charge[frame_idx])

        return g, cell, ids, np.setdiff1d(self.mmidx,sr_atom.atoms.ids)
    
    def compute_elec(self, idx = None, mmidx=None):
        if idx is not None:
            self.nc.trajectory[idx]

        results = {}
        coord_qm = self.nc.atoms.positions[self.qmidx]/CONV_UNIT_LEN

        if mmidx is None:
            coord_mm = self.nc.atoms.positions[self.mmidx]/CONV_UNIT_LEN
            mmcharges = self.nc.atoms.charges[self.mmidx]
        else:
            coord_mm = self.nc.atoms.positions[mmidx]/CONV_UNIT_LEN
            mmcharges = self.nc.atoms.charges[mmidx]

        dis_vec = coord_qm[:,None] - coord_mm[None]
        dis = np.linalg.norm(dis_vec, axis=-1)
        elec_vec = (mmcharges[None,:]/np.power(dis,3))[...,None]*dis_vec

        if self.analytical_grad:
            results['Tij1'] = elec_vec

        elec_vec = elec_vec.sum(1)
        results['T0'] = np.float32((mmcharges[None,:]/dis).sum(1))
        results['T1'] = np.float32(elec_vec)

        I = np.eye(3)
        if self.T_max_order >= 1:
            results['Tij2'] = 3*(dis_vec[:,:,:,None]*dis_vec[:,:,None,:]) - np.power(dis,2)[:,:,None,None]*I[None,None,:,:]
            results['Tij2'] = results['Tij2']/(np.power(dis,5)[...,None,None])*mmcharges[None,:,None,None]
            results['T2'] = np.float32((results['Tij2']).sum(1))

            if not (self.T_max_order ==1 and self.analytical_grad):
                results.pop('Tij2')

        if self.T_max_order >= 2: 
            results['Tij3'] = 15*dis_vec[:,:,:,None,None]*dis_vec[:,:,None,:,None]*dis_vec[:,:,None,None,:] - \
                              3*(dis**2)[:,:,None,None,None]*(dis_vec[:,:,:,None,None]*I[None,None,None,:,:] + \
                                                              dis_vec[:,:,None,:,None]*I[None,None,:,None,:] + \
                                                              dis_vec[:,:,None,None,:]*I[None,None,:,:,None])
            results['Tij3'] = (mmcharges[None,:,None,None,None]*results['Tij3']/(np.power(dis,7)[:,:,None,None,None]))
            results['T3'] = results['Tij3'].sum(1)

            if not (self.T_max_order ==2 and self.analytical_grad):
                results.pop('Tij3')
        
        if self.T_max_order >= 3:
            results['Tij4'] = 105*dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,:,None,None]*dis_vec[:,:,None,None,:,None]*dis_vec[:,:,None,None,None,:] - \
                              15*(dis**2)[:,:,None,None,None,None]*(dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,:,None,None]*I[None,None,None,None,:,:] + \
                                                                    dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,None,:,None]*I[None,None,None,:,None,:] + \
                                                                    dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,None,None,:]*I[None,None,None,:,:,None] + \
                                                                    dis_vec[:,:,None,:,None,None]*dis_vec[:,:,None,None,:,None]*I[None,None,:,None,None,:] + \
                                                                    dis_vec[:,:,None,:,None,None]*dis_vec[:,:,None,None,None,:]*I[None,None,:,None,:,None] + \
                                                                    dis_vec[:,:,None,None,:,None]*dis_vec[:,:,None,None,None,:]*I[None,None,:,:,None,None]) 
  
            results['Tij4'] = results['Tij4'] + 3*np.power(dis,4)[:,:,None,None,None,None]*(I[None,None,:,:,None,None]*I[None,None,None,None,:,:] + I[None,None,:,None,:,None]*I[None,None,None,:,None,:] + I[None,None,:,None,None,:]*I[None,None,None,:,:,None]) 
            results['Tij4'] = ((results['Tij4']/(np.power(dis,9)[:,:,None,None,None,None]))*mmcharges[None,:,None,None,None,None])
            results['T4'] = results['Tij4'].sum(1)

            if not (self.T_max_order ==3 and self.analytical_grad):
                results.pop('Tij4')

        return results

    def copy(self,nc=None,label_data=None):

        if nc is None:
            nc = self.nc

        if label_data is None:
            label_data = deepcopy(self.label_data)
            transform_unit = False
        else:
            transform_unit = True

        if 'T1' not in label_data.keys() or 'T0' not in label_data.keys():
            if_compute_elec = True 
        else:
            if_compute_elec = False
        
        new = self.__class__(datafile = nc,#.copy(),                 
                             label_file= label_data,
                             select_group=self.select_index,    
                             full_connected=self.full_connected,  
                             cutoff=self.cutoff,              
                             force_field_type=self.force_field_type,
                             add_self_loop=self.add_self_loop,   
                             if_compute_elec = if_compute_elec ,
                             T_max_order = self.T_max_order,
                             transform_unit = transform_unit
                            )
        
        if self.if_compute_elec and nc is None:
            new.T = deepcopy(self.T)
            new.if_cal = deepcopy(self.if_cal)
        return new
    
    def save(self, path = None):
        if path is None:
            label_file = 'label_data.h5'
            trj_file = 'trj.nc'
        else:
            label_file = path+'.h5'
            trj_file = path+'.dcd'

        write_data(label_file, self.label_data)

        with mda.Writer(trj_file, self.select_group.n_atoms) as w:
            for ts in self.nc.trajectory:
                w.write(self.select_group)

    def total_energy(self):
        if 'T0' in self.label_data.keys():
            totene = self.label_data['qm_energy'] + (self.label_data['qm_charge'] * self.label_data['T0'][...,None]).sum(-1).sum(-1)
        elif self.if_compute_elec and self.if_cal.sum()>0:
            totene = self.label_data['qm_energy'][self.if_cal==1] + (self.label_data['qm_charge'][self.if_cal==1] * self.T['T0'][self.if_cal==1][...,None]).sum(-1).sum(-1)        
        return totene

    def restore_dataset(self, dir):
        if self.if_cal.sum() == self.n_frames:
            for key in self.T.keys():
                self.label_data[key] = self.T[key]

        trj= mda.core.universe.Merge(self.select_group)
        filelist2 = glob.glob(dir+r"/active_epoch*.dcd")
        # trj.load_new(dataset,in_memory=True)

        label_data = defaultdict(list)
        filelist = []
        for i in range(len(filelist2)):
            labelfile = os.path.join(dir,'active_epoch='+str(i)+'.h5')
            filelist.append(os.path.join(dir,'active_epoch='+str(i)+'.dcd'))
            data = h5py.File(labelfile)

            for key in data.keys():
                label_data[key].append(data.get(key)[:])
        
        for key in label_data.keys():
            label_data[key] = th.from_numpy(np.concatenate(label_data[key]))

        trj.load_new(filelist,in_memory=True)
        coord = np.concatenate([self.nc.trajectory.coordinate_array[:,self.qmidx],
                                trj.trajectory.coordinate_array],axis=0)
        trj.load_new(coord, format=MemoryReader)
        for key in label_data.keys():
            label_data[key] = th.cat([self.label_data[key],label_data[key]])

        label_data['qm_energy'] = label_data['qm_energy']*CONV_UNIT_ENE
        label_data['qm_force'] = label_data['qm_force']*CONV_UNIT_FORCE

        dataset = self.copy(nc = trj, label_data = label_data)
        return dataset

    
class mydataset_mlp(Dataset):
    def __init__(self,
                datafile,
                label_file=None):
        self.data = np.load(datafile)['trj']
        self.label_data = np.load(label_file)
        self.cmt_forw = self.label_data['cmt_forw']
        self.cmt_back = self.label_data['cmt_back']
        self.bias = self.label_data['bias']
        self.n_frames = self.data.shape[0]

    def __getitem__(self,idx):
        if self.label_data is not None:
            return self.data[idx], self.bias[idx], self.cmt_forw[idx], self.cmt_back[idx]
        else:
            return self.data[idx]

    def __len__(self):
        return self.n_frames

def combine_trj(topfile, dir, origin_trj):
    filelist = glob.glob(dir+r"/active_epoch*.dcd")[0]
    coord = []
    for i in range(len(filelist)):
        trjfile = os.path.join(dir,'active_epoch='+str(i)+'.dcd')
        labelfile = os.path.join(dir,'active_epoch='+str(i)+'.h5')
        mda.Universe(topfile,trjfile,in_memory=True)


class graph_server():
    def __init__(self, main_rank, n_frame, n_atoms_max, n_edges_qm_max, n_atoms_qm_max):

        self.main_rank = main_rank
        self.n_atoms_max = n_atoms_max
        self.n_edges_qm_max = n_edges_qm_max
        self.n_atoms_qm_max = n_atoms_qm_max
        self.n_frame = n_frame
        self.save_time = 0
        self.done = False
        self.ndata = defaultdict(dict)
        self.edata = defaultdict(dict)
        self.edges = defaultdict(list)
        self.etype = []

        if not th.distributed.is_initialized():
            raise("to use this server, torch.distributed is must be initialized here")
        else:
            self.ddp = True

    def initialize(self,g):
        for key1 in g.ndata.keys():
            if len(g.ntypes)>1:
                for key2 in g.ndata[key1].keys():
                    if self.ddp:
                        if key2 == "mm":
                            shape1 = th.Size([self.n_frame,self.n_atoms_max])
                            shape = shape1 + g.ndata[key1][key2].shape[1:]
                            dtype = g.ndata[key1][key2].dtype
                            th.distributed.barrier()
                            self.ndata[key1][key2] =  shared_tensor(shape,dtype=dtype)
                            if self.main_rank:
                                self.ndata[key1][key2].copy_(th.zeros(shape,dtype=dtype))
                        elif key2 == "qm":
                            shape1 = th.Size([self.n_frame])
                            shape = shape1 + g.ndata[key1][key2].shape
                            dtype = g.ndata[key1][key2].dtype
                            th.distributed.barrier()
                            self.ndata[key1][key2] =  shared_tensor(shape,dtype=dtype)

                            if self.main_rank:
                                self.ndata[key1][key2].copy_(th.zeros(shape,dtype=dtype))

        if 'Hqm' in g.etypes:
            self.etype.append(('qm','Hqm','qm'))

        if 'Hqm2' in g.etypes:
            self.etype.append(('qm','Hqm2','qm'))
            self.hgraph = True
        
        if 'Hqmmm' in g.etypes:
            self.etype.append(('mm','Hqmmm','qm'))
            self.n_edges_qmmm_max = g.ndata['xyz']['qm'].shape[0]*self.n_atoms_max
            self.hgraph = True

        for key1 in g.edata.keys():
            if len(g.etypes)>1:
                for key2 in g.edata[key1].keys():
                    if len(self.etype)<len(g.etypes):
                        if key2 not in self.etype:
                            self.etype.append(key2)

                    dtype = g.edata[key1][key2].dtype
                    if key2 == ('mm','Hqmmm','qm'):
                        shape1 = th.Size([self.n_frame,self.n_edges_qmmm_max])
                        shape = shape1 + g.edata[key1][key2].shape[1:]
                    elif key2 == ('qm','Hqm','qm'):
                        if self.n_edges_qm_max > 0:
                            shape1 = th.Size([self.n_frame,self.n_edges_qm_max])
                            shape = shape1 + g.edata[key1][key2].shape[1:]
                        else:
                            shape1 = th.Size([self.n_frame])
                            shape = shape1 + g.edata[key1][key2].shape

                    th.distributed.barrier()
                    self.edata[key1][key2] = shared_tensor(shape,dtype=dtype)

                    if self.main_rank:
                        self.edata[key1][key2].copy_(th.zeros(shape,dtype=dtype))

        for key in self.etype:
            if key == ('mm','Hqmmm','qm'):
                shape = th.Size([self.n_frame,self.n_edges_qmmm_max,2])
            elif key == ('qm','Hqm','qm'):
                if self.n_edges_qm_max == 0:
                    shape1 = th.Size([self.n_frame])
                    shape = shape1 + g.edges(etype=key)[0].shape+th.Size([2])
                else:
                    shape = th.Size([self.n_frame,self.n_edges_qm_max,2])
            elif key == ('qm','Hqm2','qm'):
                n_full_edge_qm= self.n_atoms_qm_max * (self.n_atoms_qm_max - 1)
                shape = th.Size([self.n_frame,n_full_edge_qm,2])

            th.distributed.barrier()
            self.edges[key] = shared_tensor(shape,dtype=th.int64)

            if self.main_rank:
                self.edges[key].copy_(th.zeros(shape,dtype=th.int64))
        
        th.distributed.barrier()
        self.atom_num_list = shared_tensor([self.n_frame],dtype=th.int32)
        self.edge_num_list = shared_tensor([self.n_frame],dtype=th.int32)
        self.full_edge_num_list = shared_tensor([self.n_frame],dtype=th.int32)
        if self.n_edges_qm_max != 0:
            self.edge_num_list_qm = shared_tensor([self.n_frame],dtype=th.int32)
        if self.main_rank:
            self.atom_num_list.copy_(th.zeros(self.n_frame))
            self.edge_num_list.copy_(th.zeros(self.n_frame))
            self.full_edge_num_list.copy_(th.zeros(self.n_frame))

            if self.n_edges_qm_max != 0:
                self.edge_num_list_qm.copy_(th.zeros(self.n_frame))
    
    def update(self, state_dict):
        if self.main_rank:
            for key in state_dict:
                if isinstance(state_dict[key],defaultdict):
                    for key2 in state_dict[key]:
                        if key == 'ndata':
                            for key3 in state_dict[key][key2]:
                                if self.ndata[key2][key3].shape == state_dict[key][key2][key3].shape:
                                    self.ndata[key2][key3].copy_(state_dict[key][key2][key3])
                                else:
                                    raise Exception(f"shape mis-matching for {key2}")
                        elif key == 'edata':
                            for key3 in state_dict[key][key2]:
                                if self.edata[key2][key3].shape == state_dict[key][key2][key3].shape:
                                    self.edata[key2][key3].copy_(state_dict[key][key2][key3])
                                else:
                                    raise Exception(f"shape mis-matching for {key2}")
                        elif key == 'edges':
                            if self.edges[key2].shape == state_dict[key][key2].shape:
                                self.edges[key2].copy_(state_dict[key][key2])
                            else:
                                raise Exception(f"shape mis-matching for {key2}")

                elif isinstance(state_dict[key],th.Tensor):
                    if key == 'atom_num_list':
                        if self.atom_num_list.shape == state_dict[key].shape:
                            self.atom_num_list.copy_(state_dict[key])
                        else:
                            raise Exception(f"shape mis-matching for {key}")
                    elif key == 'edge_num_list':
                        if self.edge_num_list.shape == state_dict[key].shape:
                            self.edge_num_list.copy_(state_dict[key])
                        else:
                            raise Exception(f"shape mis-matching for {key}")
                    elif key == 'edge_num_list_qm':
                        if self.edge_num_list_qm.shape == state_dict[key].shape:
                            self.edge_num_list_qm.copy_(state_dict[key])
                        else:
                            raise Exception(f"shape mis-matching for {key}")
                    elif key == 'full_edge_num_list':
                        if self.full_edge_num_list.shape == state_dict[key].shape:
                            self.full_edge_num_list.copy_(state_dict[key])
                        else:
                            raise Exception(f"shape mis-matching for {key}")
        self.done = True
        self.save_time = 1
        th.distributed.barrier()

    def get_graph(self,idx):

        if self.atom_num_list[idx] != 0:
            if self.hgraph:
                eg = {}
                n_atom = self.atom_num_list[idx]
                n_full_edge_qm = self.full_edge_num_list[idx]
                n_edge = self.edge_num_list[idx]

                if self.n_edges_qm_max != 0:
                    n_edge_qm = self.edge_num_list_qm[idx]

                for key in self.edges:
                    if key == ('mm','Hqmmm','qm'):
                        eg[key] = (self.edges[key][idx,:n_edge,0],self.edges[key][idx,:n_edge,1])
                    elif key == ('qm','Hqm','qm'):
                        if self.n_edges_qm_max == 0:
                            eg[key] = (self.edges[key][idx,:,0],self.edges[key][idx,:,1])
                        else:
                            eg[key] = (self.edges[key][idx,:n_edge_qm,0],self.edges[key][idx,:n_edge_qm,1])
                    elif key == ('qm','Hqm2','qm'):
                        eg[key] = (self.edges[key][idx,:n_full_edge_qm,0],self.edges[key][idx,:n_full_edge_qm,1])
                g = dgl.heterograph(eg)

                for key1 in self.ndata.keys():
                    ndata = {}
                    for key2 in self.ndata[key1].keys():
                        if key2 == 'mm':
                            ndata[key2] = self.ndata[key1][key2][idx,:n_atom]
                        else:
                            ndata[key2] = self.ndata[key1][key2][idx]
                    g.ndata[key1] = ndata

                for key1 in self.edata.keys():
                    edata = {}
                    for key2 in self.edata[key1].keys():
                        if key == ('mm','Hqmmm','qm'):
                            edata[key2] = self.edata[key1][key2][idx,:n_edge]
                        elif key == ('qm','Hqm','qm'):
                            if self.n_edges_qm_max == 0:
                                edata[key2] = self.edata[key1][key2][idx]
                            else:
                                edata[key2] = self.edata[key1][key2][idx,:n_edge_qm]
                    g.edata[key1] = edata
            
            return g
        else:
            return None

    def append(self,g,idx):
        n_atom = g.batch_num_nodes(ntype='mm')[0]
        n_atom_qm = g.batch_num_nodes(ntype='qm')[0]
        n_edge = g.batch_num_edges(etype='Hqmmm')[0]
        if self.n_edges_qm_max > 0:
            n_edge_qm = g.batch_num_edges(etype='Hqm')[0]

            if n_edge_qm > self.n_edges_qm_max:
                raise(f"the size of truncated qm edge: {n_edge_qm} is larger than n_edges_qm_max, please set n_edges_qm_max to larger value")

        if n_atom > self.n_atoms_max or n_edge > self.n_edges_qmmm_max:
            raise(f"the size of truncated cluster: {n_atom} is larger than n_atoms_max, please set n_atoms_max to larger value")

        if n_atom_qm > self.n_atoms_qm_max:
            raise(f"the size of truncated cluster: {n_atom_qm} is larger than n_atoms_qm_max, please set n_atoms_qm_max to larger value")

        self.atom_num_list[idx].copy_(n_atom)
        self.edge_num_list[idx].copy_(n_edge)
        self.full_edge_num_list[idx].copy_(n_atom_qm*(n_atom_qm-1))

        if self.n_edges_qm_max != 0:
            self.edge_num_list_qm[idx].copy_(n_edge_qm)

        for key1 in g.ndata.keys():
            if len(g.ntypes)>1:
                for key2 in g.ndata[key1].keys():
                    if key2 == 'mm':
                        self.ndata[key1][key2][idx,:n_atom].copy_(g.ndata[key1][key2])
                    else:
                        self.ndata[key1][key2][idx].copy_(g.ndata[key1][key2])

        if 'Hqm' in g.etypes:
            self.etype.append(('qm','Hqm','qm'))
        if 'Hqm2' in g.etypes:
            self.etype.append(('qm','Hqm2','qm'))

        for key1 in g.edata.keys():
            if len(g.etypes)>1:
                for key2 in g.edata[key1].keys():

                    if key2 == ('mm', 'Hqmmm', 'qm'):
                        self.edata[key1][key2][idx,:n_edge].copy_(g.edata[key1][key2])
                    elif key2 == ('qm','Hqm','qm'):
                        if self.n_edges_qm_max == 0:
                            self.edata[key1][key2][idx].copy_(g.edata[key1][key2])
                        else:
                            self.edata[key1][key2][idx,:n_edge_qm].copy_(g.edata[key1][key2])

        for key in self.etype:
            if key == ('mm','Hqmmm','qm'):
                self.edges[key][idx,:n_edge,0].copy_(g.edges(etype=key)[0])
                self.edges[key][idx,:n_edge,1].copy_(g.edges(etype=key)[1])
            elif key == ('qm','Hqm','qm'):
                if self.n_edges_qm_max == 0 :
                    self.edges[key][idx,:,0].copy_(g.edges(etype=key)[0])
                    self.edges[key][idx,:,1].copy_(g.edges(etype=key)[1])
                else:
                    self.edges[key][idx,:n_edge_qm,0].copy_(g.edges(etype=key)[0])
                    self.edges[key][idx,:n_edge_qm,1].copy_(g.edges(etype=key)[1])
            elif key == ('qm','Hqm2','qm'):
                n_full_edge_qm = n_atom_qm * (n_atom_qm - 1)
                self.edges[key][idx,:n_full_edge_qm,0].copy_(g.edges(etype=key)[0])
                self.edges[key][idx,:n_full_edge_qm,1].copy_(g.edges(etype=key)[1])
        
        if (self.atom_num_list == 0).sum() == 0:
            self.done = True

    def save(self, filename):
        if self.main_rank: 
            state_dict = self.__dict__
            th.save(state_dict,filename)
        self.save_time = self.save_time + 1


        
        



