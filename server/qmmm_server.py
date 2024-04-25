from networkx.classes.function import neighbors
import numpy as np
import torch
import MDAnalysis.lib.nsgrid as ngrid
import dgl
from dgl.data.utils import save_graphs
import logging

class QMMMDriver:

    def __init__(self):
        logger = logging.getLogger("Driver")
        logger.info("qmmmdriver")
        self.results = None

    def calculate(self, namd_data):
        self._create_input(namd_data)
        self._run_calculation()
        results = self._parse_output()
        return results

    def _create_input(self, namd_data):
        raise NotImplementedError

    def _run_calculation(self):
        raise NotImplementedError

    def _parse_output(self):
        raise NotImplementedError

class MLDriver(QMMMDriver):

    def __init__(self, model_path, device=torch.device('cpu'),
                 fully_connected=True, cutoff=3.0):
        logger = logging.getLogger("main.Driver")
        logger.info('enter the Driver')
        super(MLDriver, self).__init__()

        self.device = device
        self.cutoff = cutoff
        self.fully_connected = fully_connected

        # Load the model
        logger.info('Loaded model...')
        self.model = self._load_model(model_path)
        self.model.eval()
        logger.info('Loaded model done')

    def _load_model(self, model_path):
        """
        Load model and activate properties for atomic charges.
        """
        model = torch.load(model_path, map_location=self.device).to(self.device)

        return model

    def _create_input(self, coor, Z, T, c, cell=None,  dtype=np.float32 ):

        n_atoms = coor.shape[0]
        if not self.fully_connected:
        ### neighbor list ###
            if cell is None:
                cell = np.array([999.0,999.0,999.0,90.0,90.0,90.0])
            else:
                if cell[:3].sum()==0:
                    cell = np.array([999.0,999.0,999.0,90.0,90.0,90.0])

            neighbor_ob=ngrid.FastNS(self.cutoff,np.float32(coor),cell)
            neighbor_list=neighbor_ob.self_search().get_pairs()
            snode = np.concatenate([neighbor_list[:,0],neighbor_list[:,1]])
            tnode = np.concatenate([neighbor_list[:,1],neighbor_list[:,0]])
            #self.g = dgl.graph((snode,tnode))
            self.cell = torch.from_numpy(cell).to(self.device)

            tnode2 = np.array([range(n_atoms)] * n_atoms)
            tnode2 = tnode2[~np.eye(tnode2.shape[0], dtype=bool)].reshape(
                                tnode2.shape[0], -1
                               )
            tnode2 = torch.from_numpy(tnode2).reshape(-1)
            snode2 = (np.repeat(np.arange(n_atoms)[:,None],n_atoms-1,axis=1)).reshape(-1)
            data_dict = {
                            ('qm', 'Hqm', 'qm'): (snode, tnode),
                            ('qm', 'Hqm2', 'qm'): (snode2, tnode2),
                        }
            self.g=dgl.heterograph(data_dict)
        else:
            ### construct graph ###
            tnode = np.array([range(n_atoms)] * n_atoms)
            tnode = tnode[~np.eye(tnode.shape[0], dtype=bool)].reshape(
                                tnode.shape[0], -1
                               )
            tnode = torch.from_numpy(tnode).reshape(-1)
            snode = (np.repeat(np.arange(n_atoms)[:,None],n_atoms-1,axis=1)).reshape(-1)
            self.cell = torch.Tensor([999.0,999.0,999.0,90.0,90.0,90.0]).to(self.device)
            self.g = dgl.graph((snode,tnode))
        
        #c = np.array([-0.21899999678134918 ,-0.11400000005960464 , -1.0, 0.11100000143051147 , 0.11100000143051147, 0.11100000143051147])
        if dtype==np.float32:
            self.g.ndata['xyz'] = torch.from_numpy(coor).float()
            self.g.ndata['z'] = torch.from_numpy(np.int32(Z)).float()
            self.g.ndata['T0'] = torch.from_numpy(T['T0']).float()
            self.g.ndata['T1'] = torch.from_numpy(T['T1']).float()
            self.g.ndata['c'] = torch.from_numpy(c).float()
            if 'T2' in T.keys():
                self.g.ndata['T2'] = torch.from_numpy(T['T2']).float()
            if 'T3' in T.keys():
                self.g.ndata['T3'] = torch.from_numpy(T['T3']).float()
            if 'T4' in T.keys():
                self.g.ndata['T4'] = torch.from_numpy(T['T4']).float()
            self.cell = self.cell.float()
        else:
            self.g.ndata['xyz'] = torch.from_numpy(coor).double()
            self.g.ndata['z'] = torch.from_numpy(np.int32(Z)).double()
            self.g.ndata['T0'] = torch.from_numpy(T['T0']).double()
            self.g.ndata['T1'] = torch.from_numpy(T['T1']).double()
            self.g.ndata['c'] = torch.from_numpy(c).double()
            if 'T2' in T.keys():
                self.g.ndata['T2'] = torch.from_numpy(T['T2']).double()
            if 'T3' in T.keys():
                self.g.ndata['T3'] = torch.from_numpy(T['T3']).double()
            if 'T4' in T.keys():
                self.g.ndata['T4'] = torch.from_numpy(T['T4']).double()
            self.cell = self.cell.double()
        #save_graphs('graph.bin',[self.g])
        self.g = self.g.to(self.device)

        #if cell is not None:
            #self.cell = torch.from_numpy(np.float32(cell)).to(self.device)


    def _run_calculation(self,dtype=np.float32):
        logger = logging.getLogger('main.calculation')
        # TODO: Activate charge derivatives!!!!
        #results = self.model(self.g,self.g.ndata['extfield'][None],cell = self.cell[None])
        results = self.model(self.g,cell = self.cell[None])
        #logger.debug(results['qm_force'].shape)
        #logger.debug(results['ext_force'].shape)
        #logger.debug(results['qm_energy'].shape)
        #logger.debug(results['ext_energy'].shape)
#        logger.debug(results['qm_charge'])
        results['qm_force'] = results['qm_force'] + results['ext_force']
        results['qm_energy'] = results['qm_energy'] + results['ext_energy']
        #logger.debug(results['qm_energy'])
        #logger.debug(results['qm_force'])
        #results_bak = self.model(self.g,cell = self.cell[None])
        #logger.info(f"{results['qm_energy'] - results_bak['qm_energy']}")
        
        
        #logger.info(results['qm_charge'].shape)
        #logger.info(self.g.ndata['extfield'].shape)
        #results['qm_force'] = results['qm_force'] + results['qm_charge'].T*self.g.ndata['extfield']

        # Convert to numpy arrays
        self.results = {}
        for prop in results:
            if dtype==np.float32:
                self.results[prop] = np.float64(results[prop].detach().cpu().numpy())
            else:
                self.results[prop] = results[prop].detach().cpu().numpy()
        return self.results

    def _parse_output(self):
        """
        Convert results to proper NAMD units
        """
        # results = {}
        # results[Properties.energy] = self.results[Properties.energy][0, 0] * au2namd_energy
        # results[Properties.forces] = self.results[Properties.forces][0] * au2namd_forces
        # results["charges"] = self.results['charges'][0]
        # return results
        return

    @staticmethod
    def _compute_charges(dipole_derivatives):
        charges = -torch.mean(torch.diagonal(dipole_derivatives, dim1=2, dim2=3), 2)
        return charges
