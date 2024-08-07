import h5py
import numpy as np
import sys


def compare(data1,data2,mm_charge=None):

    f = (abs(data1.get('qm_force_total')[:] - data2.get('qm_force_total')[:]).mean())
    e = ((abs(data1.get('energy')[:] - data2.get('energy')[:])).mean())
    mmf = (abs((data1.get('mm_force')[:] - data2.get('mm_force')[:])).mean())
    if mm_charge is not None:
        mme = (abs((data1.get('mm_force')[:] - data2.get('mm_force')[:])/mm_charge[None,:,None]).mean())
    else:
        mme = None
    return f,e,mmf,mme

if __name__ == "__main__":
    ft = []
    e = []
    mmf = []
    mme = []
    f1 = open(sys.argv[1])
    f2 = open(sys.argv[2])
    if len(sys.argv) > 4:
        mm_charge = np.loadtxt(sys.argv[3])
    else:
        mm_charge = None
    
    for (i,j) in zip(f1,f2):
        data1 = h5py.File(i.strip())
        data2 = h5py.File(j.strip())
        print(i.strip())
        print(data1.keys())
        print(data2.keys())
        f1,e1,mmf1, mme1 = compare(data1,data2,mm_charge = mm_charge)
        print(f"energy:{e1}")
        print(f"qmforce:{f1}")
        print(f"mmforce:{mmf1}")
        ft.append(f1)
        e.append(e1)
        mmf.append(mmf1)
        mme.append(mme1)
    print(np.array(ft).mean())
    print(np.array(e).mean())
    print(np.array(mmf1).mean())
    if mm_charge is not None:
        print(np.array(mme1).mean())
       

