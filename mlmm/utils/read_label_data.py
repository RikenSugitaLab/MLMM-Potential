from tokenize import Comment
import numpy as np
import os
import glob
import linecache
from collections import defaultdict
import h5py
import sys

qminfo2genesis = 0.280028557027
CONV_UNIT_ENE = 627.5095
CONV_UNIT_LEN    = 0.52917721092
CONV_UNIT_FORCE  = CONV_UNIT_ENE / CONV_UNIT_LEN

def write_data(filename,data):
    with h5py.File(filename,'w') as f:
        for i in data: 
           f[i]=data[i]

def read_qminfo(filename):

    line_id = 1
    order = 'head -n 4 '+filename+' > '+filename+'.tmpfile'
    os.system(order)
    energy = np.loadtxt(filename+'.tmpfile',comments='QM',skiprows=1)
    energy[0] = energy[0] * CONV_UNIT_ENE
    energy[1] = energy[1] * CONV_UNIT_ENE * qminfo2genesis

    line_id = 5
    qm_natoms = int(linecache.getline(filename,8).strip('\n'))

    line_id = 9
    end_line = line_id + qm_natoms - 1 
    order = "sed -n '"+str(line_id)+','+str(end_line)+"p' "+filename+' > '+filename+'.tmpfile'
    os.system(order)
    qm_charge = np.loadtxt(filename+'.tmpfile')

    line_id = end_line + 3
    end_line = line_id + qm_natoms - 1
    order = "sed -n '"+str(line_id)+','+str(end_line)+"p' "+filename+' > '+filename+'.tmpfile'
    os.system(order)
    qm_force_total = np.loadtxt(filename+'.tmpfile')
    qm_force_total = qm_force_total * CONV_UNIT_FORCE

    line_id = end_line + 3
    end_line = line_id + qm_natoms - 1
    order = "sed -n '"+str(line_id)+','+str(end_line)+"p' "+filename+' > '+filename+'.tmpfile'
    os.system(order)
    qmmm_force = np.loadtxt(filename+'.tmpfile')
    qmmm_force = qmmm_force * qminfo2genesis * CONV_UNIT_FORCE

    mm_natoms = int(linecache.getline(filename,end_line+2))
    line_id = end_line + 3
    end_line = line_id + mm_natoms - 1
    order = "sed -n '"+str(line_id)+','+str(end_line)+"p' "+filename+' > '+filename+'.tmpfile'
    os.system(order)
    mm_force = np.loadtxt(filename+'.tmpfile')
    mm_force = mm_force * CONV_UNIT_FORCE
    return energy, qm_charge, qm_force_total, qmmm_force, mm_force
    
def read_label_data(dir,outfile='label.h5py'):
    os.chdir(dir)
    file_id = glob.glob('*.qminfo')
    file_id_sort = sorted(file_id, key=lambda name: int(name[3:-7]))

    data = defaultdict(list)
    for file in file_id_sort:
        e, c, qf, qmf, mf= read_qminfo(file)
        data['energy'].append(e)
        data['qm_charge'].append(c[:,1:])
        data['qm_force_total'].append(qf[:,1:])
        data['qmmm_force'].append(qmf[:,1:])
        data['mm_force'].append(mf[:,1:])

    for key in data:
        data[key] = np.stack(data[key])
    

    label = {}
    label['extfield'] = data['qmmm_force']/data['charge'][:,:,None]
    label['qm_force'] = data['qm_force_total'] - data['qmmm_force']
    label['charge'] = data['charge']
    write_data('qminfo_'+outfile,data)
    write_data('label_'+outfile,data)
    return

if __name__=='__main__':
    read_label_data(sys.argv[1],outfile=sys.argv[2])
    #read_label_data('production.0',outfile='replica_0.h5py')
    #read_qminfo('job0.qminfo')
    #read_label_data('./')

