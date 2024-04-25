from collections import defaultdict
from multiprocessing.sharedctypes import RawArray
import pickle
from multiprocessing import synchronize
from multiprocessing.connection import Listener
from tracemalloc import get_traced_memory
import time
import torch as th
import logging
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
import dask
import dask.multiprocessing
from dask.distributed import Client
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
from multiprocessing.managers import BaseManager
class QueueManager(BaseManager): pass

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
                 num_workers = 0,
                 num_gpu = 1
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
        self.transform = None
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
        self.num_workers = num_workers
        self.num_gpu = num_gpu
        # self.authkey = os.urandom(20)
        # th.multiprocessing.set_start_method('spawn')


        self.trj_format = ['dcd','nc','filelist','npz']

        ### reading the trajectory file ###
        if isinstance(datafile, str):
            self.format = datafile.split('.')[-1]
        else:
            self.format = ''

    def load_data(self):
        self.ddp = th.distributed.is_initialized()
        
        if self.ddp:
            self.main_rank = (th.distributed.get_rank()==0)
        else:
            self.main_rank = True

        if self.format in self.trj_format and not isinstance(self.datafile, mda.Universe):

            if self.format == 'filelist':
               self.datafile = read_filelist(self.datafile)

            if self.format != 'npz':
                self.nc = mda.Universe(self.top,self.datafile,in_memory=self.in_memory,forces=True)
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
                element={'F':9,'O':8,'C':6,'N':7,'S':16,'H':1,'CL':35}
                atom_type= self.select_group.atoms.names

            self.atom_index=np.array([element[re.sub(r'[0-9]+', '', i)] for i in atom_type],dtype=np.float32)
            self.atom_mass=np.array(self.select_group.atoms.masses)

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
                       if self.sharing_data and ( self.ddp or self.num_workers>0):
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
                       if self.sharing_data and ( self.ddp or self.num_workers>0):
                            if self.ddp:
                                th.distributed.barrier()
                            self.label_data['T0'] = shared_tensor((self.n_frames,self.n_atoms_qm))
                            if self.main_rank:
                                self.label_data['T0'].copy_(th.from_numpy(f.get(key)[:]).float())
                       else:
                           self.label_data['T0'] = th.from_numpy(f.get(key)[:]).float()
                       key = 'T0'
                   else:
                       if self.sharing_data and ( self.ddp or self.num_workers>0):
                            if self.ddp:
                                th.distributed.barrier()
                            self.label_data[key] = shared_tensor(f[key].shape)
                            if self.main_rank:
                                self.label_data[key].copy_(th.from_numpy(f.get(key)[:]).float())
                       else:
                           self.label_data[key] = th.from_numpy(f.get(key)[:]).float()
                   #print(self.label_data[key].shape[0])
                   #print(self.n_frames)
                   assert self.label_data[key].shape[0]==self.n_frames, 'shape mismatch'
            else:
               raise Exception('can not process the label file ended with'+self.label_file.split('.')[-1])
        elif not isinstance(self.label_file, str):
            self.label_data = self.label_file
        else:
            self.label_data = None

        if isinstance(self.label_data, dict) and self.transform_unit:
            self.label_data['qm_energy'] = self.label_data['qm_energy']/CONV_UNIT_ENE
            self.label_data['qm_force'] = self.label_data['qm_force']/CONV_UNIT_FORCE
            if 'mm_force' in self.label_data.keys():
                if self.main_rank:
                    # self.label_data['mm_force'] = self.label_data['mm_force']/CONV_UNIT_FORCE
                    self.label_data['mm_force'].copy_(self.label_data['mm_force']/CONV_UNIT_FORCE)
            # self.label_data['T1'] = self.label_data['T1']/CONV_UNIT_FORCE

        self.get_statistic() 

        if self.num_workers > 0:
            self.barrier = mp.Barrier(self.num_workers)#*self.num_gpu)
            self.queue = mp.Queue(self.num_workers - 1)

        self.graph_server = graph_server(self.main_rank,self.n_frames,barrier=self.barrier, num_workers = self.num_workers)

        self.qmidx = np.sort(self.qmidx)
        self.mmidx = np.sort(self.mmidx)
        self.internal_idx = np.zeros(self.nc.atoms.n_atoms,dtype=np.int64)
        self.internal_idx[self.qmidx] = np.arange(self.qmidx.shape[0],dtype=np.int64) 
        self.internal_idx[self.mmidx] = np.arange(self.mmidx.shape[0],dtype=np.int64) 
        self.conn = None
        
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

    def my_barrier(self):
        curProc = mp.current_process()
        if curProc._identity[0] == 1:
            if self.main_rank:
                self.conn.send(['1'])
            else:
                synchronize = False
                while (not synchronize):
                    try:
                        self.conn.recv()
                        synchronize = True
                    except:
                        a = 2
        
    def gather_mp(self,input):
        curProc = mp.current_process()

        if curProc._identity[0] != 1:
            input_clone = input
            self.queue.put(input_clone)
        
        self.barrier.wait()

        recv_list = []
        for i in range(self.queue.qsize()):
            if curProc._identity[0] == 1:
                recv_list.append(self.queue.get())

        self.barrier.wait()
        return recv_list
    
    def __getitem__(self,idx):

        curProc = mp.current_process()

        # if self.conn is None and curProc._identity[0] == 1:
        if False: #self.conn is None and curProc._identity[0] == 1:
            address = ('localhost',6002)
            if self.main_rank:
                self.listener = mp.connection.Listener(address,authkey=b'12345')#self.authkey)
                self.conn = self.listener.accept()

            if not self.main_rank:
                no_connection = True
                while (no_connection):
                    try:
                        self.conn = mp.connection.Client(address,authkey=b'12345')#self.authkey)
                        no_connection = False
                        print("connect")
                    except:
                        print('no connection, try again')

        self.barrier.wait()

        self.nc.trajectory[idx]

        ### the frames are transformed into graph data on the fly when training ###

        # g = self.graph_server.get_graph(idx)
        g = None

        if g is None:
            g,cell,eid = self.nglist_search(idx,
                                  self.select_group,
                                  atom_index=self.atom_index,
                                  atom_mass=self.atom_mass,
                                  cutoff=self.cutoff,
                                  full_conneted=self.full_connected)

            if self.add_self_loop:
                g = g.add_self_loop()

            T = self.compute_elec()

            for key in T.keys():
                if "Tij" not in key:
                    if len(g.ntypes)<1:
                        g.ndata[key] = T[key]
                    else:
                        g.ndata[key] = {'qm': th.from_numpy(T[key]).float()}
                else:
                    g.edata[key] = {'Hqmmm':th.from_numpy(T[key][eid[0],eid[1]]).float()}
            
            # if not self.ddp:
                # self.graph_server.append(g,idx)
        else:
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

        # if self.num_workers > 0:
            # g_list_mp = self.gather_mp(g)
            # idx_list_mp = self.gather_mp(idx)
# 
            # if curProc._identity[0] == 1:
                # g_list_mp.append(g)
                # idx_list_mp.append(idx)
# 
            # if curProc._identity[0] == 1:
                # if not self.main_rank:
                    # self.conn.send(pickle.dumps(g_list_mp))
                    # print('send done')
                # else:
                    # no_message = True
                    # while (no_message):
                        # try:
                            # g_list = pickle.loads(self.conn.recv())
                            # no_message = False
                        # except:
                            # a = 2
                    # print('receive done')
            # self.my_barrier()
            # self.barrier.wait()
# 
            # if curProc._identity[0] == 1:
                # if not self.main_rank:
                    # self.conn.send(idx_list_mp)
                    # print('send done')
                # else:
                    # no_message = True
                    # while (no_message):
                        # try:
                            # idx_list = self.conn.recv()
                            # no_message = False
                        # except:
                            # a = 2
                    # print('receive done')
            # self.my_barrier()
            # self.barrier.wait()
# 
            # if curProc._identity[0] == 1 and self.main_rank:
                    # g_list.extend(g_list_mp)
                    # idx_list.extend(idx_list_mp)
            # else:
                # g_list = [1 for _ in range(self.num_gpu*self.num_workers)]
                # idx_list = [None for _ in range(self.num_gpu*self.num_workers)]
# 
            # if (self.graph_server.idx>-1).sum() != self.n_frames:
                # for arg in zip(g_list,idx_list):
                    # if arg[0] is not None:
                        # self.my_barrier()
                        # self.barrier.wait()
                        # if self.main_rank and curProc._identity[0] ==1:
                            # self.graph_server.append(*arg)
                        # else:
                            # self.graph_server.append(g,idx)
            # else:
                # if g is None:
                    # raise Exception('g is None')
        print(g.ndata['xyz']['mm'].shape[0])

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
                # cell = th.Tensor([999.0,999.0,999.0,90.0,90.0,90.0])
                # atom_groups.dimensions = cell.numpy()
                raise Exception('no periodic box information, full_connected should be true')
        ### neighbor list ###
            if atom_groups.dimensions[:3].sum()==0:
                cell = th.Tensor([999.0,999.0,999.0,90.0,90.0,90.0])
            else:
                cell = th.from_numpy(atom_groups.dimensions)

            neighbor_ob=ngrid.FastNS(cutoff,atom_groups.positions,cell)
            neighbor_list=neighbor_ob.self_search().get_pairs()
            snode = neighbor_list[:,0]
            snode = neighbor_list[:,1]
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
            sr_neighbor_ob = AtomNeighborSearch(self.nc.atoms[self.mmidx],box = cell)
            sr_atom = sr_neighbor_ob.search(atom_groups,self.sr_raidus)
            n_atoms_sr = sr_atom.n_atoms

            tnode_qm = (th.arange(n_atoms)[:,None]).repeat(1,n_atoms_sr).reshape(-1)
            snode_mm = ((th.arange(n_atoms_sr)[:,None]).repeat(1,n_atoms).T).reshape(-1)
            data_dict = {
                            ('qm', 'Hqm', 'qm'): (snode, tnode),
                            ('mm', 'Hqmmm', 'qm'): (snode_mm,tnode_qm),
                        }
            ids = (self.internal_idx[self.qmidx[tnode_qm]],self.internal_idx[sr_atom.atoms.ids[snode_mm]])
            g=dgl.heterograph(data_dict)
            g.ndata['c'] = {'mm':th.from_numpy(np.float32(sr_atom.charges))}
            g.ndata['xyz'] = {'mm': th.from_numpy(sr_atom.positions), 'qm': th.from_numpy(atom_groups.positions)}

            if hasattr(sr_atom,"force"):
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

            if hasattr(sr_atom,"force"):
                g.ndata['f'] = th.from_numpy(atom_groups.forces)
            else:
                g.ndata['f'] = self.label_data['qm_force']
            g.ndata['m'] = th.from_numpy(atom_mass)
            ids = None

            #cell = th.from_numpy(cell)

        # if charge is not None:
            # g.ndata['charge'] = th.from_numpy(charge[frame_idx])

        return g, cell, ids
    
    def compute_elec(self, idx = None):
        if idx is not None:
            self.nc.trajectory[idx]

        results = {}
        coord_qm = self.nc.atoms.positions[self.qmidx]/CONV_UNIT_LEN
        coord_mm = self.nc.atoms.positions[self.mmidx]/CONV_UNIT_LEN
        mmcharges = self.nc.atoms.charges[self.mmidx]

        dis_vec = coord_qm[:,None] - coord_mm[None]
        dis = np.linalg.norm(dis_vec, axis=-1)
        elec_vec = (mmcharges[None,:]/np.power(dis,3))[...,None]*dis_vec

        if self.analytical_grad:
            results['Tij1'] = elec_vec

        elec_vec = elec_vec.sum(1)
        results['T0'] = np.float32((mmcharges[None,:]/dis).sum(1))
        results['T1'] = np.float32(elec_vec)

        I = np.eye(3)
        if self.T_max_order > 1 or (self.T_max_order==1 and self.analytical_grad):
            results['Tij2'] = 3*(dis_vec[:,:,:,None]*dis_vec[:,:,None,:]) - np.power(dis,2)[:,:,None,None]*I[None,None,:,:]
            results['Tij2'] = results['Tij2']/(np.power(dis,5)[...,None,None])*mmcharges[None,:,None,None]

            if self.T_max_order > 1:
                results['T2'] = np.float32((results['Tij2']).sum(1))

        if self.T_max_order > 2 or (self.analytical_grad and self.T_max_order ==2):
            results['Tij3'] = 15*dis_vec[:,:,None,None,None]*dis_vec[:,:,None,:,None]*dis_vec[:,:,None,None,:] - \
                              3*(dis**2)[:,None,None,None]*(dis_vec[:,:,:,None,None]*I[None,None,None,:,:] + \
                                                          dis_vec[:,:,None,:,None]*I[None,None,:,None,:] + \
                                                          dis_vec[:,:,None,None,:]*I[None,None,:,:,None])

            results['Tij3'] = -(mmcharges[None,:,None,None,None]*results['Tij3']/(np.power(dis,5)[:,:,None,None,None]))
            results['T3'] = results['Tij3'].sum(1)
        
        if self.T_max_order == 3 and self.analytical_grad:
            results['Tij4'] = 105*dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,:,None,None]*dis_vec[:,:,None,None,:,None]*dis_vec[:,:,None,None,None,:] - \
                              15*(dis**2)[:,:,None,None,None,None]*(dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,:,None,None]*I[None,None,None,None,:,:] + \
                                                                    dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,None,:,None]*I[None,None,None,:,None,:] + \
                                                                    dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,None,None,:]*I[None,None,None,:,:,None] + \
                                                                    dis_vec[:,:,None,:,None,None]*dis_vec[:,:,None,None,:,None]*I[None,None,:,None,None,:] + \
                                                                    dis_vec[:,:,None,:,None,None]*dis_vec[:,:,None,None,None,:]*I[None,None,:,None,:,None] + \
                                                                    dis_vec[:,:,None,None,:,None]*dis_vec[:,:,None,None,None,:]*I[None,None,:,:,None,None]) 
  
            results['Tij4'] = results['Tij4'] + 3*np.power(dis,4)[:,:,None,None,None,None]*(I[None,None,:,:,None,None]*I[None,None,None,None,:,:] + I[None,None,:,None,:,None]*I[None,None,None,:,None,:] + I[None,None,:,None,None,:]*I[None,None,None,:,:,None]) 
            results['Tij4'] = ((results['Tij4']/(np.power(dis,9)[:,:,None,None,None,None]))*mmcharges[:,None,None,None,None,None])
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
                             full_connected=self.full_conneted,  
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
    def __init__(self,main_rank,n_frame,g=None,barrier = th.distributed.barrier, num_workers = 0):
        self.main_rank = main_rank
        if g is None:
            self.initialized = False
        else:
            self.initialized = True
        self.ndata = defaultdict(dict)
        self.edata = defaultdict(dict)
        self.edges = defaultdict(list)
        self.etype = []
        self.ddp = th.distributed.is_initialized()
        self.barrier = barrier
        self.num_workers = num_workers

        if self.ddp:
            th.distributed.barrier()
            self.idx = shared_tensor((n_frame,1))
            if self.main_rank:
                self.idx.copy_(th.zeros([n_frame,1])-1)
        else:
            self.idx = th.zeros((n_frame,1)) - 1
        self.internal_idx = 0

    def get_graph(self,idx2):
        idx3 = th.where(self.idx == idx2)[0]
        if idx3.shape[0] != 0:
            if self.hgraph:
                eg = {}
                for key in self.edges:
                    eg[key] = self.edges[key][idx3]
                g = dgl.heterograph(eg)

                for key1 in self.ndata.keys():
                    ndata = {}
                    for key2 in self.ndata[key1].keys():
                        ndata[key2] = self.ndata[key1][key2][idx3]
                    g.ndata[key1] = ndata

                for key1 in self.edata.keys():
                    edata = {}
                    for key2 in self.edata[key1].keys():
                        edata[key2] = self.edata[key1][key2][idx3]
                    g.edata[key1] = edata
            else:
                g = dgl.graph(self.edges[idx3])
            
            return g
        else:
            return None

    def append(self,g,idx):
        for key1 in g.ndata.keys():
            if len(g.ntypes)>1:
                for key2 in g.ndata[key1].keys():
                    if not self.initialized:
                        self.ndata[key1][key2] =  []
                        self.hgraph = True

                    if self.ddp:
                        # th.distributed.barrier()
                        self.barrier.wait()
                        self.ndata[key1][key2].append(shared_tensor(g.ndata[key1][key2].shape))
                        if self.main_rank:
                            self.ndata[key1][key2][-1].copy_(g.ndata[key1][key2])
                    else:
                        self.ndata[key1][key2].append(g.ndata[key1][key2])


        if 'Hqm' in g.etypes:
            self.etype.append(('qm','Hqm','qm'))

        for key1 in g.edata.keys():
            if len(g.etypes)>1:
                for key2 in g.edata[key1].keys():
                    if len(self.etype)<len(g.etypes):
                        if key2 not in self.etype:
                            self.etype.append(key2)
                    if not self.initialized:
                        self.edata[key1][key2] = []

                    if self.ddp:
                        # th.distributed.barrier()
                        self.barrier.wait()
                        self.edata[key1][key2].append(shared_tensor(g.edata[key1][key2].shape))

                        if self.main_rank:
                            self.edata[key1][key2][-1].copy_(g.edata[key1][key2])
                    else:
                        self.edata[key1][key2].append(g.edata[key1][key2])

        for key in self.etype:
            if self.ddp:
                # th.distributed.barrier()
                self.barrier.wait()
                self.edges[key].append((shared_tensor(g.edges(etype=key).shape)),shared_tensor(g.edges(etype=key).shape))
                if self.main_rank:
                    self.edges[key][-1][0].copy_(g.edges(etype=key)[0])
                    self.edges[key][-1][1].copy_(g.edges(etype=key)[1])
            else:
                self.edges[key].append((g.edges(etype=key)[0],g.edges(etype=key)[1]))

        if self.main_rank:
            self.idx[self.internal_idx].copy_(th.tensor([idx]))
            self.internal_idx = self.internal_idx + 1
        
        if not self.initialized:
            self.initialized = True

        



