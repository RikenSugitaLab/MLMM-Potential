import numpy as np
import torch
from collections import defaultdict
import dftbplus
import MDAnalysis as mda
import h5py

CONV_UNIT_ENE = 627.5095
CONV_UNIT_LEN    = 0.52917721092
CONV_UNIT_FORCE  = CONV_UNIT_ENE / CONV_UNIT_LEN
BOHR__AA = 0.529177249
AA__BOHR = 1 / BOHR__AA

class PotentialCalculator:
    '''

       Auxiliary class for calculating the population dependent external
       potential and its gradients. An instance of this class gets handed over
       to DFTB+ via the ctypes interface, to handle the necessary callbacks.

    '''


    def __init__(self, qmcoords, mmcoords=None, mmcharges=None, field_st = None, field = None):
        '''Initializes a PotentialCalculator object.

        Args:

            qmcoords (2darray): coordinates of QM-atoms (shape: [qmatoms, 3])
            mmcoords (2darray): coordinates of MM-atoms (shape: [mmatoms, 3])
            mmcharges (1darray): charges of MM-atoms (shape: [mmatoms, 1])

        '''

        self._qmcoords = qmcoords*AA__BOHR
        self._qmatoms = np.shape(self._qmcoords)[0]

        if mmcoords is not None:
            self._mmcoords = mmcoords*AA__BOHR
            self._mmatoms = np.shape(self._mmcoords)[0]
        else:
            self._mmcoords = mmcoords
            self._mmatoms = 0

        self._mmcharges = mmcharges
        self.field_st = field_st
        self.field = field

    def calc_extpot(self, dqatom):
        '''Calculates the current external potential using the properties of the
           MM- and QM-atoms.

        Args:

            dqatom (1darray): population difference with respect to reference
                population (usually the neutral atom). Note: population means
                electrons, so a positive number indicates electron excess

        Returns:

            extpot (1darray): updated external potential at the position of each
                QM-atom

        '''

        # Note: Some types of potential require knowledge of the
        # current atomic populations, which is provided by dqatom.

        if self.field_st is None:
            extpot = np.zeros(self._qmatoms)

            for iatqm in range(0, self._qmatoms):
                extpot[iatqm] = 0.0
                qmpos = self._qmcoords[iatqm, :]
                for iatmm in range(0, self._mmatoms):
                    mmpos = self._mmcoords[iatmm, :]
                    mmcharge = self._mmcharges[iatmm]
                    dist = distance(qmpos, mmpos)
                    extpot[iatqm] += -mmcharge / dist
        
        if self.field_st is not None:
            return -self.field_st
        else:
            return extpot


    def calc_extpotgrad(self, dqatom):
        '''Calculates the current gradients of the external potential using the
           properties of the MM- and QM-atoms.

        Args:

            dqatom (1darray): population difference with respect to reference
                population (usually the neutral atom). Note: population means
                electrons, so a positive number indicates electron excess

        Returns:

            extpotgrad (2darray): updated potential gradient at the position of
                each QM-atom

        '''

        # Note: Some types of potential require knowledge of the
        # current atomic populations, which is provided by dqatom.

        if self.field is None:
            extpotgrad = np.zeros((self._qmatoms, 3))

            for iatqm in range(0, self._qmatoms):
                qmpos = self._qmcoords[iatqm, :]
                extpotgrad[iatqm, :] = 0.0
                for iatmm in range(0, self._mmatoms):
                    mmpos = self._mmcoords[iatmm, :]
                    mmcharge = self._mmcharges[iatmm]
                    dist = distance(qmpos, mmpos)
                    dist3 = dist**3
                    extpotgrad[iatqm, 0] += -mmcharge * (mmpos[0]-qmpos[0]) / dist3
                    extpotgrad[iatqm, 1] += -mmcharge * (mmpos[1]-qmpos[1]) / dist3
                    extpotgrad[iatqm, 2] += -mmcharge * (mmpos[2]-qmpos[2]) / dist3

            return extpotgrad
        else:
            return self.field


def get_extpot(potcalc, dqatom, extpotatom):
    '''Queries the external potential.

    Args:

        potcalc (pyobject): instance of a class that provides methods for
            calculating the external potential and its gradients
        dqatom (1darray): population difference with respect to reference
            population (usually the neutral atom). Note: population means
            electrons, so a positive number indicates electron excess
        extpotatom (1darray): potential at the position of each QM-atom.
            Note: it should be the potential as felt by an electron
            (negative potential value means attraction for an electron)

    '''

    extpotatom[:] = potcalc.calc_extpot(dqatom)


def get_extpotgrad(potcalc, dqatom, extpotatomgrad):
    '''Queries the external potentials gradients.

    Args:

        potcalc (pyobject): instance of a class that provides methods for
            calculating the external potential and its gradients
        dqatom (1darray): population difference with respect to reference
            population (usually the neutral atom). Note: population means
            electrons, so a positive number indicates electron excess
        extpotatomgrad (2darray): potential gradient at the position of each
            QM-atom. Note: it should be the gradient of the potential as felt by
            an electron (negative potential value means attraction for an
            electron)

    '''

    extpotatomgrad[:, :] = potcalc.calc_extpotgrad(dqatom)


def distance(aa, bb):
    '''Auxiliary function for calculating the distance between two vectors.

    Args:

       aa (1darray): vector a
       bb (1darray): vector b

    Returns:

       dist (float): distance between vectors

    '''

    aa = np.asarray(aa)
    bb = np.asarray(bb)

    dist = np.linalg.norm(aa - bb)

    return dist


class EnergyCalculatorBase:
    def __init__(self):
        self.calculator = None
    
    def setup(self, qmcoords, mmcoords=None, mm_charges=None, field_st=None, field=None):
        self.qmcoords = qmcoords*AA__BOHR
        if mmcoords is not None:
            self.mmcoords = mmcoords*AA__BOHR
        else:
            self.mmcoords = mmcoords

        self.mm_charges = mm_charges
        self.field_st = field_st 
        self.field = field
    
    def compute(self):
        raise NotImplemented

class DftbCalculator(EnergyCalculatorBase):
    def __init__(self,lib_path, hsd_path, logfile='log.log'):
        super().__init__()
        self.lib_path = lib_path
        self.hsd_path = hsd_path
        self.logfile = logfile

    def setup(self,
              qmcoords,
              mmcoords = None,
              mm_charges = None,
              field_st = None,
              field = None
              ):
        if self.calculator is None:
            self.calculator = dftbplus.DftbPlus(libpath=self.lib_path,
                                                hsdpath=self.hsd_path,
                                                logfile=self.logfile)
        super().setup(qmcoords, mmcoords = mmcoords, mm_charges=mm_charges,field_st = field_st, field = field)
    
    def compute(self):
        if self.qmcoords.ndim  != 3:
            self.qmcoords = self.qmcoords[None]
        
        self.qmcoords = np.float64(self.qmcoords)
        if self.field is None:
            mode = 'point charge'
        else:
            mode = 'electric field'
        
        results = defaultdict(list)
        for i in range(self.qmcoords.shape[0]):

            if mode == 'point charge':
                mmcoords = self.mmcoords[i]
                field = None
                field_st = None
            else:
                mmcoords = None
                field = self.field[i]
                field_st = self.field_st[i]

            potcalc = PotentialCalculator(self.qmcoords[i], 
                                          mmcoords = mmcoords,
                                          mmcharges = self.mm_charges, 
                                          field_st = field_st,
                                          field = field)
            self.calculator.set_geometry(self.qmcoords[i],latvecs=None)
            self.calculator.register_ext_pot_generator(potcalc, get_extpot, get_extpotgrad)

            merminen = self.calculator.get_energy()
            gradients = self.calculator.get_gradients()
            grosschgs = self.calculator.get_gross_charges()

            results['qm_charge'].append(grosschgs[...,None])
            results['qm_energy'].append((merminen -  (grosschgs*self.field_st[i]).sum())*CONV_UNIT_ENE)
            results['qm_force'].append((-gradients - self.field[i]*grosschgs[...,None])*CONV_UNIT_FORCE)
        
        for key in results.keys():
            if key == 'qm_energy':
                results[key] = torch.from_numpy(np.array(results[key]))
            else:
                results[key] = torch.from_numpy(np.stack(results[key],axis=0))
        results['extfield'] = torch.from_numpy(self.field)*CONV_UNIT_FORCE
        results['extfield_st'] = torch.from_numpy(self.field_st)

        print(f"energy:{results['qm_energy'].mean()}")
        # self.calculator.close()
        return results

def compute_elec( coord_qm, coord_mm, mmcharges):
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
        
if __name__ == "__main__":

    ncfile = 'trj_al.dcd'
    label_file = 'data_al.h5'
    topfile = 'solute.pdb'

    trj = mda.Universe(topfile, ncfile, in_memory=True)
    label = h5py.File(label_file)
# 
    # field, field_st = compute_elec(trj.trajectory.coordinate_array[:2,:6],
                                #    trj.trajectory.coordinate_array[:2,6:],
                                #    trj.atoms.charges[6:])
# 
    LIB_PATH = '/home/yklei/software/anaconda3/envs/pytorch/lib/libdftbplus'
    labelor = DftbCalculator(LIB_PATH, './dftb_in.hsd')
    labelor.setup(trj.trajectory.coordinate_array,field_st = label.get('extfield_st')[:], field = label.get('extfield')[:] )
    results1 = labelor.compute()
    labelor.setup(trj.trajectory.coordinate_array,field_st = label.get('extfield_st')[:], field = label.get('extfield')[:]*CONV_UNIT_FORCE )
    results2 = labelor.compute()
    a=3




# 
        # 
# 