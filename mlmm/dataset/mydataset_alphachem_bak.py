from collections import defaultdict
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
from dask.distributed import Client
from torch.utils.data import Dataset, DataLoader
from MDAnalysis.coordinates.memory import MemoryReader
from mlmm.utils.utils import read_filelist, write_data
from copy import deepcopy
import glob
import h5py
import re
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

    # LIB_PATH = '/home/yklei/software/anaconda3/envs/pytorch/lib/libdftbplus'
    # cal = DftbCalculator(LIB_PATH,'./dftb_in.hsd')
    # cal.setup(batched_graph.ndata['xyz'].reshape(-1,6,3),field_st=labelout['extfield_st'],field=labelout['extfield'])
    # cal.setup(batched_graph.ndata['xyz'].reshape(-1,6,3).numpy(),field_st=batched_graph.ndata['extfield_st'].reshape(64,6).numpy(),field=batched_graph.ndata['extfield'].reshape(64,6,3).numpy())
    # res = cal.compute()


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
                 if_compute_elec = False, 
                 cal_quadrupole = False,
                 transform_unit = True 
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
        self.cal_quadrupole = cal_quadrupole
        self.transform_unit = transform_unit

        trj_format = ['dcd','nc','filelist','npz']

        ### reading the trajectory file ###
        if isinstance(datafile, str):
            self.format = datafile.split('.')[-1]
        else:
            self.format = ''

        if self.format in trj_format and not isinstance(self.datafile, mda.Universe):
            self.in_memory = in_memory

            if self.format == 'filelist':
               self.datafile = read_filelist(self.datafile)

            if self.format != 'npz':
                self.nc = mda.Universe(top,self.datafile,in_memory=in_memory)
                self.n_frames = self.nc.trajectory.n_frames
                if self.nc.trajectory.velocity_array is not None:
                    self.read_vel=True
            elif self.format == 'npz':
                data = np.load(self.datafile)
                if top != '':
                    self.nc = mda.Universe(top,data['R'],format=MemoryReader)
                else:
                    n_atoms = data['R'].shape[1]
                    self.nc = mda.Universe.empty(n_atoms,trajectory=True,velocities=True)
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

        self.select_group = self.nc.select_atoms(select_group)
        self.qmidx = self.select_group.ix_array
        self.mmidx = np.setdiff1d(np.arange(self.nc.atoms.n_atoms), self.qmidx)
        if top != '' or isinstance(self.datafile, mda.Universe):
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

        self.full_conneted=full_connected

        if label_file != '' and isinstance(label_file,str):
            if label_file.split('.')[-1]=='npz':
               label_data = np.load(label_file)
               self.label_data = {}
               for key in list(label_data):
                   if key in ['z','R','v']:
                       continue
                   self.label_data[key] = th.from_numpy(label_data[key])
                   assert self.label_data[key].shape[0]==self.n_frames, 'shape mismatch'
            elif label_file.split('.')[-1]=='txt':
               self.label_data =th.from_numpy(np.loadtxt(label_file))
               assert self.label_data.shape[0]==self.n_frames, 'shape mismatch'
            elif label_file.split('.')[-1]=='h5py' or label_file.split('.')[-1]== 'h5':
                f = h5py.File(label_file,'r')
                self.label_data = {}
                for key in f.keys(): 
                   if key in ['z','R','v']:
                       continue
                #    self.label_data[key] = th.from_numpy(f.get(key)[:][1:]).float()
                   self.label_data[key] = th.from_numpy(f.get(key)[:]).float()
                   #print(self.label_data[key].shape[0])
                   #print(self.n_frames)
                   assert self.label_data[key].shape[0]==self.n_frames, 'shape mismatch'
            else:
               raise Exception('can not process the label file ended with'+label_file.split('.')[-1])
        elif not isinstance(label_file, str):
            self.label_data = label_file
        else:
            self.label_data = None

        if isinstance(self.label_data, dict) and self.transform_unit:
            self.label_data['extfield'] = self.label_data['extfield']/CONV_UNIT_FORCE
            self.label_data['qm_energy'] = self.label_data['qm_energy']/CONV_UNIT_ENE
            self.label_data['qm_force'] = self.label_data['qm_force']/CONV_UNIT_FORCE
        self.get_statistic() 

        if self.if_compute_elec:
            n_atoms = self.select_group.n_atoms
            self.extfield = th.zeros([self.n_frames, n_atoms, 3]).float()
            self.extfield_st = th.zeros([self.n_frames, n_atoms]).float()
            self.quadrupole = th.zeros([self.n_frames, n_atoms,3,3]).float()
            self.need_cal = np.zeros([self.n_frames])
        
        #super(mydataset,self).__init__(name='dataset_name')
        super(mydataset,self).__init__()
                
    def get_statistic(self):
        self.mean_e = self.label_data['qm_energy'].mean()
        self.std_e = self.label_data['qm_energy'].std()

    def __getitem__(self,idx):
        ### the frames are transformed into graph data on the fly when training ###
        g,cell= self.nglist_search(idx,
                              self.select_group,
                              atom_index=self.atom_index,
                              atom_mass=self.atom_mass,
                              cutoff=self.cutoff,
                              full_conneted=self.full_conneted)

        if self.add_self_loop:
            g = g.add_self_loop()

        #if self.if_compute_elec and 'extfield_st' not in self.label_data.keys():
        if self.if_compute_elec: # and 'extfield_st' not in self.label_data.keys():
            if self.need_cal[idx] == 0:
                elec_vec, elec_st, quadrupole = self.compute_elec(idx)
                self.extfield[idx] = th.from_numpy(elec_vec)
                self.extfield_st[idx] = th.from_numpy(elec_st)
                g.ndata['extfield_st'] = th.from_numpy(elec_st)
                #test
                g.ndata['extfield'] = th.from_numpy(elec_st)
                #end
                if self.cal_quadrupole:
                    g.ndata['quadrupole'] = th.from_numpy(quadrupole)
                self.need_cal[idx] = 1
            else:
                if self.cal_quadrupole:
                    g.ndata['quadrupole'] = self.quadrupole[idx]
                g.ndata['extfield_st'] = self.extfield_st[idx]

        if self.label_data is not None:
            label = {}
            keyname = []
            for key in self.label_data:
                if key not in ['R','v','extfield','extfield_st']:

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

                if key=="extfield":
                    g.ndata['extfield'] = self.label_data[key][idx]

                if key=="extfield_st":
                    g.ndata['extfield_st'] = self.label_data[key][idx]

            return g, keyname, tuple(label.values()), cell
        else:
            return g, cell 

    def __len__(self):
        return self.n_frames
    
    def nglist_search(self,frame_idx,atom_groups,atom_index,cutoff,atom_mass=None,full_conneted=False):
        ### transform the frames into graph data using mdanalysis ###
        n_atoms = atom_groups.n_atoms
        atom_groups.universe.trajectory[frame_idx]
        coordinates = atom_groups.positions

        if self.read_vel: 
            velocities = atom_groups.velocities

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
            neighbor_ob=ngrid.FastNS(cutoff,coordinates,cell)
            neighbor_list=neighbor_ob.self_search().get_pairs()
            g=dgl.graph((th.IntTensor(neighbor_list[:,0]),th.IntTensor(neighbor_list[:,1])),num_nodes=n_atoms)
        else:
            ### construct graph ###
            tnode = np.array([range(n_atoms)] * n_atoms)
            tnode = tnode[~np.eye(tnode.shape[0], dtype=bool)].reshape(
                                tnode.shape[0], -1
                               )
            tnode = th.from_numpy(tnode).reshape(-1)
            snode = (np.repeat(np.arange(n_atoms)[:,None],n_atoms-1,axis=1)).reshape(-1)
            cell = th.Tensor([999.0,999.0,999.0,90.0,90.0,90.0])
            #cell = np.array([999.0,999.0,999.0,90.0,90.0,90.0])
            g = dgl.graph((snode,tnode))

        ### node feature ###
        g.ndata['xyz'] = th.from_numpy(coordinates)
        if self.read_vel:
            g.ndata['v'] = th.from_numpy(velocities)
        g.ndata['z'] = th.from_numpy(np.int32(atom_index))
        g.ndata['m'] = th.from_numpy(atom_mass)
        #cell = th.from_numpy(cell)

        # if charge is not None:
            # g.ndata['charge'] = th.from_numpy(charge[frame_idx])

        return g, cell
    
    def compute_elec(self, idx):
        coord_qm = self.nc.trajectory.coordinate_array[idx,self.qmidx]
        coord_mm = self.nc.trajectory.coordinate_array[idx,self.mmidx]
        mmcharges = self.nc.atoms.charges[self.mmidx]

        dis_vec = coord_qm[:,None] - coord_mm[None,]
        dis = np.linalg.norm(dis_vec, axis=-1)
        elec_vec = (mmcharges[None,:]/np.power(dis,3))[...,None]*dis_vec
        elec_vec = elec_vec.sum(1)*(CONV_UNIT_LEN**2)
        elec_st = (mmcharges[None,:]/dis).sum(1)*(CONV_UNIT_LEN)

        if self.cal_quadrupole:
            dis_vec_norm = dis_vec/dis[...,None]
            quadrupole = (dis_vec_norm[...,None]@dis_vec_norm[:,:,None])*mmcharges[None,:,None,None]
            quadrupole = quadrupole/(np.power(dis,3)[...,None,None])
            quadrupole = (quadrupole*(CONV_UNIT_LEN**3)).sum(1)
        else:
            quadrupole = None
        return np.float32(elec_vec), np.float32(elec_st), np.float32(quadrupole)

    def copy(self,nc=None,label_data=None):

        if nc is None:
            nc = self.nc

        if label_data is None:
            label_data = deepcopy(self.label_data)
            transform_unit = False
        else:
            transform_unit = True

        if 'extfield' not in label_data.keys() or 'extfield_st' not in label_data.keys():
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
                            #  transform_unit = transform_unit
                         )
        
        if self.if_compute_elec and nc is None:
            new.extfield = deepcopy(self.extfield)
            new.extfield_st = deepcopy(self.extfield_st)
            new.need_cal = deepcopy(self.need_cal)
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
        if 'extfield_st' in self.label_data.keys():
            totene = self.label_data['qm_energy'] + (self.label_data['qm_charge'] * self.label_data['extfield_st'][...,None]).sum(-1).sum(-1)
        elif self.if_compute_elec and self.need_cal.sum()>0:
            totene = self.label_data['qm_energy'][self.need_cal==1] + (self.label_data['qm_charge'][self.need_cal==1] * self.extfield_st[self.need_cal==1][...,None]).sum(-1).sum(-1)        
        return totene

    def restore_dataset(self, dir):
        if self.need_cal.sum() == self.n_frames:
            self.label_data['extfield_st'] = self.extfield_st
            

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
        if self.need_cal.sum() == self.n_frames:
            self.label_data['extfield_st'] = self.extfield_st
            coord = np.concatenate([self.nc.trajectory.coordinate_array[:,self.qmidx],
                                    trj.trajectory.coordinate_array],axis=0)
            trj.load_new(coord, format=MemoryReader)
            for key in label_data.keys():
                label_data[key] = th.cat([self.label_data[key],label_data[key]])

        label_data['extfield'] = label_data['extfield']*CONV_UNIT_FORCE
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


        

