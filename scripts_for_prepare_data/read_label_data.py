from tokenize import Comment
from multiprocessing import cpu_count
from MDAnalysis.coordinates.memory import MemoryReader
import argparse
import numpy as np
from pmda.parallel import ParallelAnalysisBase
import os
import glob
import linecache
from collections import defaultdict
import h5py
import sys
import dask
import dask.multiprocessing
from dask.distributed import Client
# from dask_jobqueue import PBSCluster
import MDAnalysis as mda
import json

#cluster = PBSCluster(cores=1,memory='30G')
dask.config.set(scheduler='processes')
os.environ["HDF5_USE_FILE_LOCKING"] = 'FALSE'

qminfo2genesis = 0.280028557027
CONV_UNIT_ENE = 627.5095
CONV_UNIT_LEN    = 0.52917721092
CONV_UNIT_FORCE  = CONV_UNIT_ENE / CONV_UNIT_LEN

def compute_ref_model(T,M):
    energy_ext = (T['T0']*M['M0']).sum(-1)
    force_ext = (T['T1']*M['M0'][...,None])
    if 'T1' in T.keys() and 'M1' in M.keys() and 'T2' in T.keys(): 
        energy_ext = energy_ext - (T['T1']*M['M1']).sum(-1).sum(-1)
        force_ext = force_ext - (T['T2']*M['M1'][...,None]).sum(1)
    if 'T2' in T.keys() and 'M2' in M.keys() and 'T3' in T.keys(): 
        energy_ext = energy_ext + (T['T2']*M['M2']).sum(-1).sum(-1).sum(-1)/2
        force_ext = force_ext + (T['T3']*M['M2'][...,None]).sum(1).sum(1)/2
    if 'T3' in T.keys() and 'M3' in M.keys() and 'T4' in T.keys(): 
        energy_ext = energy_ext - (T['T3']*M['M3']).sum(-1).sum(-1).sum(-1).sum(-1)/6
        force_ext = force_ext - (T['T4']*M['M3'][...,None]).sum(1).sum(1).sum(1)/6
    return energy_ext, force_ext

class compute_elec_parallel(ParallelAnalysisBase):
    def __init__(self, atomgroup, qmidx, mmidx, max_order):
        self.qmidx = qmidx
        self.mmidx = mmidx
        self.max_order = max_order
        self._ag = atomgroup
        super(compute_elec_parallel,self).__init__(atomgroup.universe, 
                                                   [self._ag])
    
    def _single_frame(self,ts, agroups):
        agroups = agroups[0]
        coord_qm = agroups.positions[self.qmidx]
        coord_mm = agroups.positions[self.mmidx]
        mmcharges = agroups.atoms.charges[self.mmidx]
        results = compute_elec(coord_qm,coord_mm,mmcharges,max_order=self.max_order)
        return results
    
    def _conclude(self):
        self._results = (self._results.reshape([-1])).tolist()
        self.results = defaultdict(list)
        if isinstance(self._results[0],dict):
            for res in self._results:
                for key in res.keys():
                    self.results[key].append(res[key])
        else:
            for res1 in self._results:
                for res2 in res1.tolist():
                    for key in res2.keys():
                        self.results[key].append(res2[key])

        for key in self.results.keys():
            self.results[key] = np.stack(self.results[key],axis=0)

def compute_elec(coor_qm, coor_mm, mmcharges, max_order = 1):
    results = {}

    coor_qm = coor_qm/CONV_UNIT_LEN
    coor_mm = coor_mm/CONV_UNIT_LEN
    dis_vec = coor_qm[:,None] - coor_mm[None]

    dis = np.linalg.norm(dis_vec, axis=-1)
    elec_vec = (mmcharges[None,:]/np.power(dis,3))[...,None]*dis_vec

    out={}
    results['Tij1'] = elec_vec

    elec_vec = elec_vec.sum(1)

    out['T0'] = np.float32((mmcharges[None,:]/dis).sum(1))
    out['T1'] = np.float32(elec_vec)

    I = np.eye(3)
    if max_order > 0 :
        results['Tij2'] = 3*(dis_vec[:,:,:,None]*dis_vec[:,:,None,:]) - np.power(dis,2)[:,:,None,None]*I[None,None,:,:]
        results['Tij2'] = results['Tij2']/(np.power(dis,5)[...,None,None])*mmcharges[None,:,None,None]

        out['T2'] = np.float32((results['Tij2']).sum(1))

    if max_order > 1:
        results['Tij3'] = 15*dis_vec[:,:,:,None,None]*dis_vec[:,:,None,:,None]*dis_vec[:,:,None,None,:] - \
                          3*(dis**2)[:,:,None,None,None]*(dis_vec[:,:,:,None,None]*I[None,None,None,:,:] + \
                                                          dis_vec[:,:,None,:,None]*I[None,None,:,None,:] + \
                                                          dis_vec[:,:,None,None,:]*I[None,None,:,:,None])

        results['Tij3'] = (mmcharges[None,:,None,None,None]*results['Tij3']/(np.power(dis,7)[:,:,None,None,None]))
        out['T3'] = results['Tij3'].sum(1)

    if max_order > 2:
        results['Tij4'] = 105*dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,:,None,None]*dis_vec[:,:,None,None,:,None]*dis_vec[:,:,None,None,None,:] - \
                          15*(dis**2)[:,:,None,None,None,None]*(dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,:,None,None]*I[None,None,None,None,:,:] + \
                                                                dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,None,:,None]*I[None,None,None,:,None,:] + \
                                                                dis_vec[:,:,:,None,None,None]*dis_vec[:,:,None,None,None,:]*I[None,None,None,:,:,None] + \
                                                                dis_vec[:,:,None,:,None,None]*dis_vec[:,:,None,None,:,None]*I[None,None,:,None,None,:] + \
                                                                dis_vec[:,:,None,:,None,None]*dis_vec[:,:,None,None,None,:]*I[None,None,:,None,:,None] + \
                                                                dis_vec[:,:,None,None,:,None]*dis_vec[:,:,None,None,None,:]*I[None,None,:,:,None,None]) 

        results['Tij4'] = results['Tij4'] + 3*np.power(dis,4)[:,:,None,None,None,None]*(I[None,None,:,:,None,None]*I[None,None,None,None,:,:] + I[None,None,:,None,:,None]*I[None,None,None,:,None,:] + I[None,None,:,None,None,:]*I[None,None,None,:,:,None]) 
        results['Tij4'] = ((results['Tij4']/(np.power(dis,9)[:,:,None,None,None,None]))*mmcharges[None,:,None,None,None,None])
        out['T4'] = results['Tij4'].sum(1)
    return out

def read_json(filename):
    f = open(filename,'r')
    data_json = json.load(f)
    prop = defaultdict(list)
    for atom_obj in data_json['multipole_moments']:
        for key in atom_obj.keys():
            if key == 'charge':
                prop['M0'].append(atom_obj[key])
            if key == 'mbis_dipole': 
                prop['M1'].append(atom_obj[key])
            elif key == 'mbis_qudrapole':
                prop['M2'].append(construct_M2(atom_obj[key]))
            elif key == 'mbis_octupole':
                prop['M3'].append(construct_M3(atom_obj[key]))
    
    out = {}
    for key in prop.keys():
        if key == 'M0':
            out[key] = np.array(prop[key])
        else:
            out[key] = np.stack(prop[key])
    return tuple(out.values())

def construct_M2(quadrapole):
    qpole = np.zeros((3,3))
    idx = 0
    for i in range(3):
        for j in range(i,3):
            qpole[i,j] = quadrapole[idx]
            qpole[j,i] = quadrapole[idx]
            idx += 1
    return qpole

def construct_M3(octupole):
    opole = np.zeros((3,3,3))
    idx = 0
    for i in range(3):
        for j in range(i,3):
            for k in range(j,3):
                opole[i,j,k] = octupole[idx]
                opole[i,k,j] = octupole[idx]
                opole[k,i,j] = octupole[idx]
                opole[k,j,i] = octupole[idx]
                opole[j,i,k] = octupole[idx]
                opole[j,k,i] = octupole[idx]
                idx += 1
    return opole

def write_data(filename,data):
    with h5py.File(filename,'w') as f:
        for i in data: 
           f[i]=data[i]

def compute_elec2( coord_qm, coord_mm, mmcharges):
        dis_vec = coord_qm[:,:,None] - coord_mm[:,None,]
        dis = np.linalg.norm(dis_vec, axis=-1)
        if len(mmcharges.shape) == 1:
            elec_vec = (mmcharges[None,None,:]/np.power(dis,3))[...,None]*dis_vec
        else:
            elec_vec = (mmcharges[:,None,:]/np.power(dis,3))[...,None]*dis_vec
        elec_vec = elec_vec.sum(2)*(CONV_UNIT_LEN**2)
        if len(mmcharges.shape) == 1:
            elec_st = (mmcharges[None,None,:]/dis).sum(2)*(CONV_UNIT_LEN)
        else:
            elec_st = (mmcharges[:,None,:]/dis).sum(2)*(CONV_UNIT_LEN)
        return np.float32(elec_vec), np.float32(elec_st)

def read_qminfo(filename):

    print("file:"+filename)
    # print("dir:"+os.getcwd())

    order = 'head -n 2 '+filename+' > '+filename+'.tmpfile'
    os.system(order)
    energy = np.loadtxt(filename+'.tmpfile',comments='QM',skiprows=1)
    energy = energy * CONV_UNIT_ENE
    line_id = 3

    end_line = line_id + 1
    order = "sed -n '"+str(end_line)+','+str(end_line)+"p' "+filename+' > '+filename+'.tmpfile'
    os.system(order)
    qm_dipole = np.loadtxt(filename+'.tmpfile')

    qm_natoms = int(linecache.getline(filename,6).strip('\n'))
    line_id = 7

    end_line = line_id + qm_natoms - 1 
    order = "sed -n '"+str(line_id)+','+str(end_line)+"p' "+filename+' > '+filename+'.tmpfile'
    os.system(order)
    qm_charge = np.loadtxt(filename+'.tmpfile')

    title = linecache.getline(filename,end_line+1).strip('\n')
    if title == "Link atom charge":
        link_atoms = int(linecache.getline(filename,end_line+2).strip('\n'))
        line_id = end_line + 3
        end_line = line_id + link_atoms - 1
        order = "sed -n '"+str(line_id)+','+str(end_line)+"p' "+filename+' > '+filename+'.tmpfile'
        os.system(order)
        link_charge = np.loadtxt(filename+'.tmpfile')
        link_atom = True
    else:
        link_atom = False
        link_charge = np.array([0])

    line_id = end_line + 3
    end_line = line_id + qm_natoms - 1
    order = "sed -n '"+str(line_id)+','+str(end_line)+"p' "+filename+' > '+filename+'.tmpfile'
    os.system(order)
    qm_force_total = np.loadtxt(filename+'.tmpfile')
    qm_force_total[:,1:4] = qm_force_total[:,1:4] * CONV_UNIT_FORCE

    if link_atom:
        line_id = end_line + 3
        end_line = line_id + link_atoms - 1
        order = "sed -n '"+str(line_id)+','+str(end_line)+"p' "+filename+' > '+filename+'.tmpfile'
        os.system(order)
        link_force = np.loadtxt(filename+'.tmpfile')
        link_force[:,2:5] = link_force[:,2:5] * CONV_UNIT_FORCE
    else:
        link_force = np.array([0])

    mm_natoms = int(linecache.getline(filename,end_line+2))
    line_id = end_line + 3
    end_line = line_id + mm_natoms - 1
    order = "sed -n '"+str(line_id)+','+str(end_line)+"p' "+filename+' > '+filename+'.tmpfile'
    os.system(order)
    mm_force = np.loadtxt(filename+'.tmpfile')
    mm_force[:,1:4] = mm_force[:,1:4] * CONV_UNIT_FORCE

    if os.path.exists(filename+'.json'):
        multipole = read_json(filename+'.json')
    else:
        multipole = np.array([0])
    return energy, qm_charge, qm_force_total, mm_force, qm_dipole, multipole, link_charge, link_force
    
def read_label_data(dirlist,outfile='label.h5py',dcdfile=None,topfile=None,interval=1,drop_first=True, ref_model = -1, store_nc = False,linkbond=1.0,log_file='log'):

    oring_dir = os.getcwd()
    data = defaultdict(list)

    if dcdfile is None:
       if ref_model >= 0:
          raise Exception("when dcdfile is None, ref_model should be less than 0")
    else:
       if dcdfile.split('.')[-1]=="filelist":
          inputfile = read_filelist(dcdfile)
       else:
          inputfile = dcdfile
       trj = mda.Universe(topfile,inputfile,in_memory=True) #,forces=True)
       natoms = trj.atoms.n_atoms 
       print(f"frames:{trj.trajectory.n_frames}")
    
    sys.stdout = open(log_file,'w')
    n_jobs=cpu_count()
    print(f"cpu core used: {n_jobs}")
    print(f"dirlist:{dirlist}")

    if not isinstance(dirlist,list):
        dirlist = [dirlist]  
    print('cpu core number:'+str(n_jobs))
    did = 0
    for dir in dirlist:
      print("enter the dir:"+dir)
      os.chdir(dir)
      print(os.getcwd())
      file_id = glob.glob('*.qminfo')
      file_id_sort = sorted(file_id, key=lambda name: int(name[3:-7]))
      os.chdir(oring_dir)

      num = len(file_id_sort)
      if did == 0:
        #cluster.scale(jobs=n_jobs)
        #client = Client(cluster)
        client = Client(n_workers=n_jobs)
        did = did + 1
      fid = 1
      joblist = []
      for file in file_id_sort:
          if int(file[3:-7])%interval != 0: 
             continue

          if drop_first and int(file[3:-7])==0:
             continue
          print(str(fid)+'/'+str(num)+' :'+os.path.join(dir,file))
          joblist.append(dask.delayed(read_qminfo)(os.path.join(dir,file)))
        #   joblist.append(read_qminfo(os.path.join(dir,file)))
        #   joblist.append(read_qminfo(file))
          sys.stdout.flush()
          fid = fid + 1

    result = dask.compute(joblist)
    client.close()
    print('read qminfo done')
    sys.stdout.flush()

    e, c, qf, mf, dip, multipole, link_charge, link_force = list(zip(*result[0]))
    e = np.stack(e)
    c = np.stack(c)
    qf = np.stack(qf)
    dip = np.stack(dip)
    mf = np.stack(mf)
    link_charge = np.stack(link_charge)
    link_force = np.stack(link_force)

    data['energy'] = e
    data['qm_charge'] = c[:,:,1:]
    data['qm_force_total'] = qf[:,:,1:]
    data['mm_force'] = mf[:,:,1:]
    data['dipole'] = dip

    #if link_charge.sum() != 0 :
    #    if not store_nc:
    #        raise Exception('when there is link atoms, store_nc should be set to true')
    if link_charge.sum() != 0 :
        data['link_charge'] = link_charge[:,:,2:]
        data['link_force'] = link_force[:,:,2:]

    label = {}
    multipole = list(zip(*multipole))
    for i in range(len(multipole)):
        if len(multipole[i][0].shape) == 1:
            label['M0'] = np.stack(multipole[i])
        elif len(multipole[i][0].shape) == 2:
            label['M1'] = np.stack(multipole[i])
        elif len(multipole[i][0].shape) == 3:
            label['M2'] = np.stack(multipole[i])
        elif len(multipole[i][0].shape) == 4:
            label['M3'] = np.stack(multipole[i])

    qmidx = np.int32(qf[0,:,0]) - 1
    mmidx = np.int32(mf[0,:,0]) - 1
    if store_nc:
        force = np.zeros_like(trj.trajectory.coordinate_array,dtype=np.float64)
        force[:,qmidx] = data['qm_force_total']
        force[:,mmidx] = data['mm_force']
        trj.trajectory.force_array = force

    if link_charge.sum() != 0:
        data['qm_charge'] =  np.concatenate([data['qm_charge'],link_charge[:,:,2:]],axis=1)
        #qmidx = np.concatenate([qmidx, link_charge[:,:,0]])
        if store_nc:
            num_link_atoms = link_charge.shape[1]
            link_atom_ob = mda.Universe.empty(num_link_atoms,trajectory=True,velocities=True,forces=True)
            link_atom_ob.add_TopologyAttr('names')
            link_atom_ob.add_TopologyAttr('charges')
            link_atom_ob.add_TopologyAttr('masses')
            link_atom_ob.atoms.names = np.repeat(np.array(['LA']), [num_link_atoms])
            link_atom_ob.atoms.charges = np.zeros([num_link_atoms],dtype=np.float64)
            link_atom_ob.atoms.masses = np.zeros([num_link_atoms],dtype=np.float64)
        #idx_exclude_atom = np.setdiff1d(np.arange(natoms),np.concatenate([qmidx,mmidx]))
 
            posa = trj.trajectory.coordinate_array[:,np.int32(link_charge[0,:,0] - 1)]
            posb = trj.trajectory.coordinate_array[:,np.int32(link_charge[0,:,1] - 1)]
            rab = np.linalg.norm(posa - posb, axis = -1)
            vab = posb - posa
            pos_link = vab * (linkbond/rab)[...,None] + posa
 
            link_atom_ob.load_new(pos_link,format=MemoryReader)
            link_atom_ob.trajectory.force_array = link_force[:,:,2:]
 
            trj2 = mda.core.universe.Merge(trj.atoms[np.concatenate([qmidx,mmidx])],link_atom_ob.atoms)
            force1 = trj.trajectory.force_array[:,np.concatenate([qmidx,mmidx])]
            force2 = link_atom_ob.trajectory.force_array
            coord1 = trj.trajectory.coordinate_array[:,np.concatenate([qmidx,mmidx])]
            coord2 = link_atom_ob.trajectory.coordinate_array
            trj2.load_new(np.concatenate([coord1,coord2],axis=1),format=MemoryReader)
            trj2.trajectory.force_array = np.concatenate([force1,force2],axis=1)
            print(f"n_frame:{trj2.trajectory.n_frames}")

    if store_nc:
        print('save the nc file with force')
        sys.stdout.flush()
        if link_force.sum() != 0:
            trj2.atoms.write('modified.pdbqt')
            w = mda.Writer(dcdfile[:-3]+'nc', trj2.atoms.n_atoms, forces=True)
            select_group = trj2.select_atoms('all')
            k = 0
            for ts in trj2.trajectory:
               print(k)
               sys.stdout.flush()
               w.write(select_group)
               k = k + 1
            trj2 = mda.Universe('modified.pdbqt',dcdfile[:-3]+'nc')
        else:
            w = mda.Writer(dcdfile[:-3]+'nc', trj.atoms.n_atoms, forces=True)
            select_group = trj.select_atoms('all')
            for ts in trj.trajectory:
               w.write(select_group)
        w.close()

    if ref_model >= 0:
        print('compute ref model')
        sys.stdout.flush()
        if link_force.sum() != 0:
            natoms2 = trj2.atoms.n_atoms
            qmidx2 = np.concatenate([np.arange(qmidx.shape[0]),np.arange(natoms2 - num_link_atoms, natoms2)])
            mmidx2 = np.arange(qmidx.shape[0],natoms2-num_link_atoms)
            
            #extfield,extfield_st = compute_elec(np.concatenate([trj.trajectory.coordinate_array[:,qmidx],pos_link],axis=0),trj.trajectory.coordinate_array[:,mmidx],trj.atoms.charges[mmidx])
            calc = compute_elec_parallel(trj2.atoms,qmidx=qmidx2,mmidx=mmidx2,max_order=ref_model)
        else:
            calc = compute_elec_parallel(trj.atoms,qmidx=qmidx,mmidx=mmidx,max_order=ref_model)
            #extfield,extfield_st = compute_elec(trj.trajectory.coordinate_array[:,qmidx],trj.trajectory.coordinate_array[:,mmidx],trj.atoms.charges[mmidx])
        calc.run(n_jobs=n_jobs,n_blocks=n_jobs)
        print('compute ref model done')
        sys.stdout.flush()
        M = {}
        if 'M0' in label.keys() or ref_model==0:
            M['M0'] = data['qm_charge'][:,:,0]
        else:
            M = multipole
        ref_ene, ref_force = compute_ref_model(calc.results,M)
        # extfield,extfield_st = compute_elec2(trj.trajectory.coordinate_array[:,qmidx],trj.trajectory.coordinate_array[:,mmidx],trj.atoms.charges[mmidx])
        # print(abs(calc.results['T0'] - extfield_st).mean())
        # print(abs(calc.results['T1'] - extfield).mean())
        # qmmm_energy = (data['qm_charge'].squeeze()*extfield_st).sum(-1)*CONV_UNIT_ENE
        # data['ref_force2'] = (data['qm_charge'].squeeze()[:,:,None]*extfield)*CONV_UNIT_FORCE
        # data['ref_energy2'] = qmmm_energy
        data['ref_energy'] = ref_ene*CONV_UNIT_ENE
        data['ref_force'] = ref_force*CONV_UNIT_FORCE
        label['ref_energy'] = data['ref_energy']
        label['ref_force'] = data['ref_force']
        # print(abs(data['ref_force2'] - data['ref_force']).mean())
        # print(abs(data['ref_energy2'] - data['ref_energy']).mean())

    print('save the h5 file')
    sys.stdout.flush()
    label['qm_charge'] = data['qm_charge']
    label['qm_dipole'] = data['dipole']
    label['qm_energy'] = data['energy']
    write_data('qminfo_'+outfile,data)
    write_data('label_'+outfile,label)
    return

def read_filelist(filename):
    with open(filename,'r') as f:
        pos_file=f.readlines()
  
        for i in range(len(pos_file)):
           pos_file[i]=pos_file[i].strip('\n')
    return pos_file

if __name__=='__main__':
    parser = argparse.ArgumentParser()

    ### params ###
    parser.add_argument("--dir",type=str, help="directory containting *.qminfo file output from genesis")
    parser.add_argument("--trj",type=str, help="trajectory file")
    parser.add_argument("--top",type=str, help="psf file")
    parser.add_argument("--outfile",type=str, help="the name output h5 file")
    parser.add_argument("--drop_first",action='store_true',default=False,help="whether job0.qminfo will be processed" )
    parser.add_argument("--store_nc",action='store_true',default=False, help="whether store coordinate and force into nc file")
    parser.add_argument("--ref_model",type=int, default=0, help="whether calculate energy and force from coloumb model, 0 is yes, -1 is no")
    parser.add_argument("--interval",type=int,default=1, help="interval for processing qminfo file, when (frame idx)%inteval==0 the corresponding qminfo file will processed")
    parser.add_argument("--linkbond",type=float,default=1.0, help="bond lenght of link bond")
    parser.add_argument("--logfile",type=str,default='log')

    args = parser.parse_args()
    print(args.dir)

    if args.dir.split('.')[-1]=="filelist":
      dirlist = read_filelist(args.dir)
      print("filelist:")
      print(dirlist)
      print(args.ref_model)
      read_label_data(dirlist,outfile=args.outfile,drop_first=args.drop_first,interval=args.interval,topfile=args.top,dcdfile=args.trj,ref_model = args.ref_model, store_nc = args.store_nc,linkbond=args.linkbond,log_file=args.logfile)
    else:
      print("dir:"+args.dir)
      read_label_data(args.dir,outfile=args.outfile,drop_first=args.drop_first,interval=args.interval,topfile=args.top,dcdfile=args.trj, ref_model = args.ref_model, store_nc = args.store_nc,linkbond=args.linkbond)
