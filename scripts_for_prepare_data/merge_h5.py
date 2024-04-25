import h5py
import sys
from collections import defaultdict
import numpy as np

def write_data(filename,data):
    with h5py.File(filename,'w') as f:
        for i in data:
           f[i]=data[i]

def merge(filelist,outfile='out.h5'):
    data_all = defaultdict(list)
    idx=0
    for i in filelist:
        file = h5py.File(i,'r')
        print(i)

        for k in file.keys():
            data_all[k].append(file.get(k)[:])

    out = {}
    for k in file.keys():
        out[k] = np.concatenate(data_all[k],axis=0)

    write_data(outfile,out)

def read_filelist(filename):
    with open(filename,'r') as f:
        pos_file=f.readlines()

        for i in range(len(pos_file)):
           pos_file[i]=pos_file[i].strip('\n')
    return pos_file

if __name__ == "__main__":
    flist = read_filelist(sys.argv[1])
    #flist = read_filelist('label.filelist')
    merge(flist,outfile=sys.argv[2])
    #merge(flist,outfile='label.h5')
    


