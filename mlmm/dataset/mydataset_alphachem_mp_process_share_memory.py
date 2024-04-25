import chunk
from collections import defaultdict
import multiprocessing
from unittest import result
import torch as th
import logging
import os
import dgl
from dgl import save_graphs, load_graphs
import netCDF4 as nc
from dgl.data.utils import makedirs, save_info, load_info
from dgl.data import DGLDataset
import MDAnalysis as mda
import numpy as np
import MDAnalysis.lib.nsgrid as ngrid
from multiprocessing import cpu_count
from multiprocessing import Pool
from functools import partial
import dask
import dask.multiprocessing
import warnings
from dask.distributed import Client
# from dask_jobqueue import PBSCluster
from torch.utils.data import Dataset, DataLoader
from MDAnalysis.coordinates.memory import MemoryReader
from mlmm.utils.utils import read_filelist, write_data
from copy import deepcopy
import glob
import h5py
import re
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
from MDAnalysis.lib.distances import capped_distance
import ctypes
import multiprocessing as mp
# from MDAnalysis import pmdm
# from mlmm.baal.get_label import DftbCalculator
from MDAnalysis.lib.util import unique_int_1d
import dask.array as da
from memory_profiler import profile
import psutil

qminfo2genesis = 0.280028557027
CONV_UNIT_ENE = 627.5095
CONV_UNIT_LEN    = 0.52917721092
CONV_UNIT_FORCE  = CONV_UNIT_ENE / CONV_UNIT_LEN
os.environ["HDF5_USE_FILE_LOCKING"] = 'FALSE'

def get_current_memory_gb():
# 获取当前进程内存占用。
    pid = os.getpid()
    p = psutil.Process(pid)
    info = p.memory_full_info()
    return info.uss / 1024. / 1024. / 1024.

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

def init_worker(coor_,
                vel_,
                force_,
                charges_,
                dimension_,
                qmidx_,
                mmidx_,
                atom_index_,
                atom_mass_,
                cutoff_,
                sr_raidus_,
                full_conneted_
                ):

    global coor 
    coor = coor_
    global vel
    vel = vel_
    global forces 
    forces = force_
    global charge 
    charge = charges_
    global dimension
    dimension = dimension_
    global qmidx
    qmidx = qmidx_
    global mmidx
    mmidx = mmidx_
    global atom_index
    atom_index = atom_index_
    global atom_mass
    atom_mass = atom_mass_
    global cutoff
    cutoff = cutoff_
    global sr_raidus
    sr_raidus = sr_raidus_
    global full_conneted
    full_conneted = full_conneted_
    


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
                 T_max_order = 1,
                 analytical_grad = True,
                 preprocess = True,
                 n_workers = 6,
                 graph_file = 'test.bin', 
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
        self.top = top
        self.label_file = label_file
        self.read_vel = False
        self.in_memory = in_memory
        self.full_conneted = full_connected
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
        self.preprocess = preprocess
        self.n_workers = n_workers
        self.graph_file = graph_file

        trj_format = ['dcd','nc','filelist','npz']

        ### reading the trajectory file ###
        if isinstance(self.datafile, str):
            self.format = self.datafile.split('.')[-1]
        else:
            self.format = ''

        if self.format in trj_format and not isinstance(self.datafile, mda.Universe):
            self.in_memory = self.in_memory

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

        self.select_group = self.nc.select_atoms(self.select_index)
        self.qmidx = self.select_group.ix_array
        self.mmidx = np.setdiff1d(np.arange(self.nc.atoms.n_atoms), self.qmidx)
        self.n_atoms_mm = self.mmidx.shape[0]
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
                       self.label_data['T1'] = th.from_numpy(f.get(key)[:]).float()
                       key = 'T1'
                   elif key == "extfield_st":
                       self.label_data['T0'] = th.from_numpy(f.get(key)[:]).float()
                       key = 'T0'
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
        
        f.close()

        if isinstance(self.label_data, dict) and self.transform_unit:
            self.label_data['qm_energy'] = self.label_data['qm_energy']/CONV_UNIT_ENE
            self.label_data['qm_force'] = self.label_data['qm_force']/CONV_UNIT_FORCE
            if 'mm_force' in self.label_data.keys():
                self.label_data['mm_force'] = self.label_data['mm_force']/CONV_UNIT_FORCE
            self.label_data['T1'] = self.label_data['T1']/CONV_UNIT_FORCE
            
        if "mm_force" in self.label_data.keys():
            # self.nc.trajectory.force_array = da.from_array(np.zeros([self.n_frames,self.nc.trajectory.n_atoms,3]),chunks=(1,self.nc.trajectory.n_atoms,3))
            self.nc.trajectory.force_array = np.zeros([self.n_frames,self.nc.trajectory.n_atoms,3])
            self.nc.trajectory.force_array[:,self.qmidx] = self.label_data['qm_force'].numpy()
            self.nc.trajectory.force_array[:,self.mmidx] = self.label_data['mm_force'].numpy()
        # self.nc.trajectory.force_array.compute()

        if 'mm_force' in self.label_data.keys():
            self.label_data.pop('mm_force')

        self.get_statistic() 

        self.graph = []
        self.if_cal = np.zeros([self.n_frames])
        self.idx = []

        self.qmidx = np.sort(self.qmidx)
        self.mmidx = np.sort(self.mmidx)
        self.internal_idx = np.zeros(self.nc.atoms.n_atoms,dtype=np.int64)
        self.internal_idx[self.qmidx] = np.arange(self.qmidx.shape[0],dtype=np.int64) 
        self.internal_idx[self.mmidx] = np.arange(self.mmidx.shape[0],dtype=np.int64) 
        self.prepare_graphlist(n_workers=self.n_workers)

        super(mydataset,self).__init__()
    
    @profile
    def prepare_graphlist(self,n_workers=None):
        # from dask.diagnostics import ResourceProfiler
        # rprof = ResourceProfiler(dt=0.5)

        # if n_workers == 1:
            # dask.config.set(scheduler='synchronous')
        # else:
            # dask.config.set(scheduler='processes')

        # client = Client(n_workers=n_workers)
        print("start")
        if n_workers is None:
            n_workers = multiprocessing.cpu_count()

        if self.preprocess and self.in_memory:
            ctx = multiprocessing.get_context('spawn')

            pool  = multiprocessing.pool.Pool(n_workers,initializer=init_worker,initargs=(self.nc.trajectory.coordinate_array,
                                                    None,
                                                    self.nc.trajectory.force_array,
                                                    self.nc.atoms.charges,
                                                    None,
                                                    self.qmidx,
                                                    self.mmidx,
                                                    self.atom_index,
                                                    self.atom_mass,
                                                    self.cutoff,
                                                    self.sr_raidus,
                                                    self.full_conneted
                                                    ),context = ctx) 
            # pool.map_async(nglist_search2,range(0,self.n_frames,1))
            async_results = [pool.apply_async(nglist_search2, args=[i]) for i in range(self.n_frames)]
            results = [ar.get() for ar in async_results]
            pool.close()
            pool.join()

            # with ResourceProfiler(dt=0.25) as rprof:
                # results = dask.compute(joblist)
# 
            # joblist = []
            # for i in range(self.n_frames):
#                # print(f"frame: {i}")
##                joblist.append(nglist_search(self.nc.trajectory.coordinate_array[i],
                # coor = dask.delayed(self.nc.trajectory.coordinate_array[i])
                # force = dask.delayed(self.nc.trajectory.force_array[i])
                # joblist.append(dask.delayed(nglist_search)(coor,#self.nc.trajectory.coordinate_array[i],
                                                        #    None,
                                                        #    force, #self.nc.trajectory.force_array[i],
                                                        #    self.nc.atoms.charges,
                                                        #    self.qmidx,
                                                        #    self.mmidx,
                                                        #    i,
                                                        #    atom_index=self.atom_index,
                                                        #    atom_mass=self.atom_mass,
                                                        #    cutoff=self.cutoff,
                                                        #    sr_raidus = self.sr_raidus,
                                                        #    full_conneted=self.full_conneted))
            # with ResourceProfiler(dt=0.25) as rprof:
                # results = dask.compute(joblist)
# 
            if self.graph_file is not None:
                keys = list(self.label_data.keys())
                for key in keys:
                    if key not in ['qm_energy','qm_charge']:
                        self.label_data.pop(key)
                dgl.save_graphs(self.graph_file,results[0],labels=self.label_data)

            return {'graph':results[0],'label':self.label_data,'cell':self.nc.trajectory.dimensions_array}

                
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

    def __getitem__(self,idx):
        if not self.preprocess:
            self.nc.trajectory[idx]
        ### the frames are transformed into graph data on the fly when training ###

        if self.if_cal[idx] == 0 or self.preprocess:
            g,cell,eid = self.nglist_search(idx,
                                  self.select_group,
                                  atom_index=self.atom_index,
                                  atom_mass=self.atom_mass,
                                  cutoff=self.cutoff,
                                  full_conneted=self.full_conneted)

            if not self.preprocess:
                T = self.compute_elec()
            else:
                T = self.compute_elec(idx=idx)

            if not self.preprocess:
                self.idx.append(idx)
            for key in T.keys():
                if "Tij" not in key:
                    if len(g.ntypes)<1:
                        g.ndata[key] = T[key]
                    else:
                        g.ndata[key] = {'qm': th.from_numpy(T[key]).float()}
                else:
                    g.edata[key] = {'Hqmmm':th.from_numpy(T[key][eid[0],eid[1]]).float()}
        else:
            idx2 = np.where(np.array(self.idx) == idx)[0][0]

            if self.nc.atoms.dimensions is None:
                cell = th.Tensor([999.0,999.0,999.0,90.0,90.0,90.0])
            else:
                cell = self.nc.atoms.dimensions

            g = self.graph[idx2]

        if self.label_data is not None and not self.preprocess:
            label = {}
            keyname = []
            for key in self.label_data:
                if key not in ['R','v','T0','T1','T2','T3','T4','mm_force']:

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

        if self.if_cal[idx] == 0 and not self.preprocess:
            self.if_cal[idx] = 1
            self.graph.append(g)

        if self.label_data is not None and not self.preprocess:
            return g, keyname, tuple(label.values()), cell
        elif self.preprocess:
            return g
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
                if not self.preprocess:
                    cell = th.from_numpy(atom_groups.dimensions)
                    neighbor_ob=ngrid.FastNS(cutoff,atom_groups.positions,cell)
                else:
                    cell = th.from_numpy(self.nc.trajectory.dimension_array[frame_idx])
                    neighbor_ob=ngrid.FastNS(cutoff,self.nc.trajectory.coordinate_array[frame_idx],cell)

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
            if not self.preprocess:
                sr_neighbor_ob = AtomNeighborSearch(self.nc.atoms[self.mmidx],box = cell)
                sr_atom = sr_neighbor_ob.search(atom_groups,self.sr_raidus)
                n_atoms_sr = sr_atom.n_atoms
            else:
                sr_atom_list =  capped_distance(self.nc.trajectory.coordinate_array[frame_idx,self.qmidx], \
                                                self.nc.trajectory.coordinate_array[frame_idx,self.mmidx], \
                                                self.sr_raidus, \
                                                box = cell, \
                                                return_distances=False)
                if sr_atom_list.size > 0:
                    sr_atom_list = self.mmidx[unique_int_1d(np.asarray(sr_atom_list[:, 1], dtype=np.intp))]
                n_atoms_sr = sr_atom_list.shape[0]

            tnode_qm = (th.arange(n_atoms)[:,None]).repeat(1,n_atoms_sr).reshape(-1)
            snode_mm = ((th.arange(n_atoms_sr)[:,None]).repeat(1,n_atoms).T).reshape(-1)
            data_dict = {
                            ('qm', 'Hqm', 'qm'): (snode, tnode),
                            ('mm', 'Hqmmm', 'qm'): (snode_mm,tnode_qm),
                        }
            
            if not self.preprocess:
                ids = (self.internal_idx[self.qmidx[tnode_qm]],self.internal_idx[sr_atom.atoms.ids[snode_mm]])
            else:
                ids = (self.internal_idx[self.qmidx[tnode_qm]],self.internal_idx[sr_atom_list[snode_mm]])
            g=dgl.heterograph(data_dict)

            if not self.preprocess:
                g.ndata['c'] = {'mm':th.from_numpy(np.float32(sr_atom.charges))}
                g.ndata['xyz'] = {'mm': th.from_numpy(sr_atom.positions), 'qm': th.from_numpy(atom_groups.positions)}
                g.ndata['f'] = {'mm': th.from_numpy(sr_atom.forces), 'qm': th.from_numpy(atom_groups.forces)}
                if self.read_vel:
                    g.ndata['v'] = {'mm': th.from_numpy(sr_atom.velocities), 'qm': th.from_numpy(atom_groups.velocities)}
            else:
                g.ndata['c'] = {'mm':th.from_numpy(np.float32(self.nc.atoms.charges[sr_atom_list]))}
                g.ndata['xyz'] = {'mm': th.from_numpy(self.nc.trajectory.coordinate_array[frame_idx,sr_atom_list]),  
                                  'qm': th.from_numpy(self.nc.trajectory.coordinate_array[frame_idx,self.qmidx])}
                g.ndata['f'] = {'mm': th.from_numpy(self.nc.trajectory.force_array[frame_idx,sr_atom_list]), 
                                'qm': th.from_numpy(self.nc.trajectory.force_array[frame_idx,self.qmidx])}
                if self.read_vel:
                    g.ndata['v'] = {'mm': th.from_numpy(self.nc.trajectory.velocity_array[frame_idx,sr_atom_list]),
                                    'qm': th.from_numpy(self.nc.trajectory.velocity_array[frame_idx,self.qmidx])}
            g.ndata['z'] = {'qm': th.from_numpy(np.int32(atom_index))}
            g.ndata['m'] = {'qm': th.from_numpy(atom_mass)}
        else:
            g=dgl.graph((snode,tnode))
            ### node feature ###
            if not self.preprocess:
                g.ndata['xyz'] = th.from_numpy(atom_groups.coordinates)
                if self.read_vel:
                    g.ndata['v'] = th.from_numpy(atom_groups.velocities)
                g.ndata['f'] = th.from_numpy(atom_groups.forces)
            else:
                g.ndata['xyz'] = th.from_numpy(self.nc.trajectory.coordinate_array[frame_idx,self.qmidx])
                if self.read_vel:
                    g.ndata['v'] = th.from_numpy(self.nc.trajectory.velocity_array[frame_idx,self.qmidx])
                g.ndata['f'] = th.from_numpy(self.nc.trajectory.force_array[frame_idx,self.qmidx])

            g.ndata['m'] = th.from_numpy(atom_mass)
            g.ndata['z'] = th.from_numpy(np.int32(atom_index))
            ids = None

        if self.add_self_loop:
            g = g.add_self_loop()

            #cell = th.from_numpy(cell)

        # if charge is not None:
            # g.ndata['charge'] = th.from_numpy(charge[frame_idx])

        return g, cell, ids
    
    def compute_elec(self, idx = None):

        results = {}
        if idx is None:
            coord_qm = self.nc.atoms.positions[self.qmidx]/CONV_UNIT_LEN
            coord_mm = self.nc.atoms.positions[self.mmidx]/CONV_UNIT_LEN
        else:
            coord_qm = self.nc.trajectory.coordinate_array[idx,self.qmidx]/CONV_UNIT_LEN
            coord_mm = self.nc.trajectory.coordinate_array[idx,self.mmidx]/CONV_UNIT_LEN
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


def nglist_search(atom_groups,
                  velocity,
                  force,
                  charge,
                  qmidx,
                  mmidx,
                  idx,
                  atom_index,
                  cutoff,
                  sr_raidus = 0,
                  dimension = None,
                  atom_mass=None,
                  full_conneted=False):
    ### transform the frames into graph data using mdanalysis ###
    print(idx)
    n_atoms = qmidx.shape[0]
    qmidx = np.sort(qmidx)
    mmidx = np.sort(mmidx)
    internal_idx = np.zeros(atom_groups.shape[0],dtype=np.int64)
    internal_idx[qmidx] = np.arange(qmidx.shape[0],dtype=np.int64) 
    internal_idx[mmidx] = np.arange(mmidx.shape[0],dtype=np.int64) 

    if velocity is not None:
        read_vel = True
    else:
        read_vel = False

    if not full_conneted:
        if dimension is None:
            # cell = th.Tensor([999.0,999.0,999.0,90.0,90.0,90.0])
            # atom_groups.dimensions = cell.numpy()
            raise Exception('no periodic box information, full_connected should be true')
    ### neighbor list ###
        cell = th.from_numpy(dimension)
        neighbor_ob=ngrid.FastNS(cutoff,atom_groups,cell)

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


    if sr_raidus > 0.0:
        sr_atom_list =  capped_distance(atom_groups[qmidx], \
                                        atom_groups[mmidx], \
                                        sr_raidus, \
                                        box = cell, \
                                        return_distances=False)
        if sr_atom_list.size > 0:
            sr_atom_list = mmidx[unique_int_1d(np.asarray(sr_atom_list[:, 1], dtype=np.intp))]
        n_atoms_sr = sr_atom_list.shape[0]

        tnode_qm = (th.arange(n_atoms)[:,None]).repeat(1,n_atoms_sr).reshape(-1)
        snode_mm = ((th.arange(n_atoms_sr)[:,None]).repeat(1,n_atoms).T).reshape(-1)
        data_dict = {
                        ('qm', 'Hqm', 'qm'): (snode, tnode),
                        ('mm', 'Hqmmm', 'qm'): (snode_mm,tnode_qm),
                    }
        
        ids = (internal_idx[qmidx[tnode_qm]],internal_idx[sr_atom_list[snode_mm]])
        g=dgl.heterograph(data_dict)

        g.ndata['c'] = {'mm':th.from_numpy(np.float32(charge[sr_atom_list]))}
        g.ndata['xyz'] = {'mm': th.from_numpy(atom_groups[sr_atom_list]),  
                          'qm': th.from_numpy(atom_groups[qmidx])}
        g.ndata['f'] = {'mm': th.from_numpy(np.array(force[sr_atom_list])), 
                        'qm': th.from_numpy(np.array(force[qmidx]))}
        if read_vel:
            g.ndata['v'] = {'mm': th.from_numpy(velocity[sr_atom_list]),
                            'qm': th.from_numpy(velocity[qmidx])}
        g.ndata['z'] = {'qm': th.from_numpy(np.int32(atom_index))}
        g.ndata['m'] = {'qm': th.from_numpy(atom_mass)}
    else:
        g=dgl.graph((snode,tnode))
        ### node feature ###
        g.ndata['xyz'] = th.from_numpy(atom_groups)
        if read_vel:
            g.ndata['v'] = th.from_numpy(velocity)
        g.ndata['f'] = th.from_numpy(force)
        g.ndata['m'] = th.from_numpy(atom_mass)
        g.ndata['z'] = th.from_numpy(np.int32(atom_index))
        ids = None
    
    T = compute_elec(atom_groups, charge, qmidx, mmidx)

    for key in T.keys():
        if "Tij" not in key:
            if len(g.ntypes)<1:
                g.ndata[key] = T[key]
            else:
                g.ndata[key] = {'qm': th.from_numpy(T[key]).float()}
        else:
            g.edata[key] = {'Hqmmm':th.from_numpy(T[key][ids[0],ids[1]]).float()}

    # if self.add_self_loop:
        # g = g.add_self_loop()

    return g #, cell, ids

def compute_elec(positions,charges,qmidx,mmidx,
                 analytical_grad=True,T_max_order=2):

    results = {}
    coord_qm = positions[qmidx]/CONV_UNIT_LEN
    coord_mm = positions[mmidx]/CONV_UNIT_LEN
    mmcharges = charges[mmidx]

    dis_vec = coord_qm[:,None] - coord_mm[None]
    dis = np.linalg.norm(dis_vec, axis=-1)
    elec_vec = (mmcharges[None,:]/np.power(dis,3))[...,None]*dis_vec

    if analytical_grad:
        results['Tij1'] = elec_vec

    elec_vec = elec_vec.sum(1)
    results['T0'] = np.float32((mmcharges[None,:]/dis).sum(1))
    results['T1'] = np.float32(elec_vec)

    I = np.eye(3)
    if T_max_order > 1 or (T_max_order==1 and analytical_grad):
        results['Tij2'] = 3*(dis_vec[:,:,:,None]*dis_vec[:,:,None,:]) - np.power(dis,2)[:,:,None,None]*I[None,None,:,:]
        results['Tij2'] = results['Tij2']/(np.power(dis,5)[...,None,None])*mmcharges[None,:,None,None]

        if T_max_order > 1:
            results['T2'] = np.float32((results['Tij2']).sum(1))

    if T_max_order > 2 or (analytical_grad and T_max_order ==2):
        results['Tij3'] = 15*dis_vec[:,:,:,None,None]*dis_vec[:,:,None,:,None]*dis_vec[:,:,None,None,:] - \
                          3*(dis**2)[:,:,None,None,None]*(dis_vec[:,:,:,None,None]*I[None,None,None,:,:] + \
                                                      dis_vec[:,:,None,:,None]*I[None,None,:,None,:] + \
                                                      dis_vec[:,:,None,None,:]*I[None,None,:,:,None])

        results['Tij3'] = -(mmcharges[None,:,None,None,None]*results['Tij3']/(np.power(dis,5)[:,:,None,None,None]))
        results['T3'] = results['Tij3'].sum(1)
    
    if T_max_order == 3 and analytical_grad:
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

def nglist_search2(idx):
    ### transform the frames into graph data using mdanalysis ###
    atom_groups = coor[idx]
    force = forces[idx]

    if vel is not None:
        velocity = vel[idx]
    else:
        velocity = None

    n_atoms = qmidx.shape[0]
    # qmidx = np.sort(qmidx)
    # mmidx = np.sort(mmidx)
    internal_idx = np.zeros(atom_groups.shape[0],dtype=np.int64)
    internal_idx[qmidx] = np.arange(qmidx.shape[0],dtype=np.int64) 
    internal_idx[mmidx] = np.arange(mmidx.shape[0],dtype=np.int64) 

    if velocity is not None:
        read_vel = True
    else:
        read_vel = False

    if not full_conneted:
        if dimension is None:
            # cell = th.Tensor([999.0,999.0,999.0,90.0,90.0,90.0])
            # atom_groups.dimensions = cell.numpy()
            raise Exception('no periodic box information, full_connected should be true')
    ### neighbor list ###
        cell = th.from_numpy(dimension)
        neighbor_ob=ngrid.FastNS(cutoff,atom_groups,cell)

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


    if sr_raidus > 0.0:
        sr_atom_list =  capped_distance(atom_groups[qmidx], \
                                        atom_groups[mmidx], \
                                        sr_raidus, \
                                        box = cell, \
                                        return_distances=False)
        if sr_atom_list.size > 0:
            sr_atom_list = mmidx[unique_int_1d(np.asarray(sr_atom_list[:, 1], dtype=np.intp))]
        n_atoms_sr = sr_atom_list.shape[0]
 
        tnode_qm = (th.arange(n_atoms)[:,None]).repeat(1,n_atoms_sr).reshape(-1)
        snode_mm = ((th.arange(n_atoms_sr)[:,None]).repeat(1,n_atoms).T).reshape(-1)
        data_dict = {
                        ('qm', 'Hqm', 'qm'): (snode, tnode),
                        ('mm', 'Hqmmm', 'qm'): (snode_mm,tnode_qm),
                    }
         
        ids = (internal_idx[qmidx[tnode_qm]],internal_idx[sr_atom_list[snode_mm]])
        g=dgl.heterograph(data_dict)
        g.ndata['c'] = {'mm':th.from_numpy(np.float32(charge[sr_atom_list]))}
        g.ndata['xyz'] = {'mm': th.from_numpy(atom_groups[sr_atom_list]),  
                          'qm': th.from_numpy(atom_groups[qmidx])}
        g.ndata['f'] = {'mm': th.from_numpy(force[sr_atom_list]), 
                        'qm': th.from_numpy(force[qmidx])}
        if read_vel:
            g.ndata['v'] = {'mm': th.from_numpy(velocity[sr_atom_list]),
                            'qm': th.from_numpy(velocity[qmidx])}
        g.ndata['z'] = {'qm': th.from_numpy(np.int32(atom_index))}
        g.ndata['m'] = {'qm': th.from_numpy(atom_mass)}
    else:
        g=dgl.graph((snode,tnode))
        ### node feature ###
        g.ndata['xyz'] = th.from_numpy(atom_groups)
        if read_vel:
            g.ndata['v'] = th.from_numpy(velocity)
        g.ndata['f'] = th.from_numpy(force)
        g.ndata['m'] = th.from_numpy(atom_mass)
        g.ndata['z'] = th.from_numpy(np.int32(atom_index))
        ids = None
    
    T = compute_elec(atom_groups, charge, qmidx, mmidx)

    for key in T.keys():
        if "Tij" not in key:
            if len(g.ntypes)<1:
                g.ndata[key] = T[key]
            else:
                g.ndata[key] = {'qm': th.from_numpy(T[key]).float()}
        else:
            g.edata[key] = {'Hqmmm':th.from_numpy(T[key][ids[0],ids[1]]).float()}
            a = 1

    # # if self.add_self_loop:
    #     # g = g.add_self_loop()
    print(idx)
    return g #, cell, ids
