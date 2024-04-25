import torch
import dgl
import pytorch_lightning as pl
import os
from torch.autograd import grad
from dgl.nn.pytorch.glob import SumPooling
from torch.autograd import Variable
import numpy as np
from mlmm.utils.utils import write_data, compute_interacton_tensor,compute_mm_analytical_grad
from mlmm.opt.diffmod import DiffMod
from mlmm.dataset.mydataset_alphachem import CONV_UNIT_ENE, CONV_UNIT_FORCE, CONV_UNIT_LEN

__all__ = ['LitMLMM']
class LitMLMM(pl.LightningModule):
    def __init__(self, rep_model: torch.nn.Module, out_model: torch.nn.Module,
                 wf: float,
                 wq: float,
                 wl2: float,
                 wg: float = 0.0,
                 wr: float = 0.0,
                 reference_model = True,
                 train_total_energy = False,
                 total_charge = -1.0,
                 mbis: bool = True):
        super().__init__()
        self.rep_model = rep_model
        self.out_model = out_model
        self.mbis = mbis
        self.train_total_energy = train_total_energy
        self.reference_model = reference_model
        self.total_charge = total_charge

        if out_model.name == "gp":
            self.mll = None
        self.pooling1 = dgl.nn.pytorch.glob.SumPooling()
        self.pooling2 = dgl.nn.pytorch.glob.AvgPooling()

        mask2 = torch.FloatTensor([[1,0,0],[1,1,0],[1,1,1]])
        mask3 = torch.zeros((3,3,3))

        for i in range(3):
            for j in range(i,3):
                for k in range(j,3):
                    mask3[i,j,k] = 1.0

        self.register_buffer('std_e',torch.FloatTensor([1.0]))
        self.register_buffer('mean_e',torch.FloatTensor([0.0]))
        if self.mbis:
            self.register_buffer('mask_M2',mask2)
            self.register_buffer('mask_M3',mask3)
        else:
            self.mask_M2 = mask2
            self.mask_M3 = mask3
        self.wf = wf
        self.wq = wq
        self.wl2 = wl2
        self.wg = wg
        self.wr = wr

    def forward(self, g, cell=None):
        results = {}
        g.ndata['xyz'].requires_grad_()

        T_keys = []
        max_order = 0
        for key in g.ndata.keys():
            if key in ['xyz','T0','T1','T2','T3']:
                g.ndata[key].requires_grad_()
                if key != 'xyz':
                    if max_order < int(key[-1]):
                        max_order = int(key[-1])
                    T_keys.append(key)
        
        if not ('T3' in g.ndata.keys() and 'T4' in g.ndata.keys()):
            T_keys.remove('T'+str(max_order))

        n_atom = g.batch_num_nodes()[0]
        n_frame = g.batch_size

        if self.out_model.name == 'mlp':
            feat, feat_v, graph_rep = self.rep_model(g,cell=cell)
            multipole, energy= self.out_model(g, feat, feat_v, graph_rep = graph_rep)
            charge = multipole['M0'].reshape(-1,n_atom)
            charge = charge + (self.total_charge - charge.sum(-1))[:,None]/n_atom
            energy = (energy * self.std_e) + self.mean_e

            if self.reference_model:
                energy_ext =  (charge*g.ndata['T0'].reshape(-1,n_atom)).sum(-1)[:,None]
                force_ext = (g.ndata['T1'].reshape(-1,n_atom,3)*(charge.squeeze())[...,None])

                if 'M1' in multipole.keys() and 'T1' in g.ndata.keys() and 'T2' in g.ndata.keys():
                    multipole['M1'] = multipole['M1'].reshape(n_frame,n_atom,3)
                    energy_ext =  energy_ext - ((multipole['M1']*g.ndata['T1'].reshape(n_frame,n_atom,3)).sum(-1).sum(-1))[:,None]
                    force_ext = force_ext - (g.ndata['T2'].reshape(-1,n_atom,3,3)*multipole['M1'][...,None]).sum(2)

                if 'M2' in multipole.keys() and 'T2' in g.ndata.keys() and 'T3' in g.ndata.keys():
                    multipole['M2'] = multipole['M2'].reshape(n_frame,n_atom,3,3)
                    energy_ext =  energy_ext + 0.5*((multipole['M2']*g.ndata['T2'].reshape(n_frame,n_atom,3,3)).sum(-1).sum(-1).sum(-1))[:,None]
                    force_ext = force_ext + (g.ndata['T3'].reshape(-1,n_atom,3,3,3)*multipole['M2'][...,None]).sum(2).sum(2)/2

                if 'M3' in multipole.keys() and 'T3' in g.ndata.keys() and 'T4' in g.ndata.keys():
                    print('multipole_M3')
                    multipole['M3'] = multipole['M3'].reshape(n_frame,n_atom,3,3,3)
                    energy_ext =  energy_ext - ((multipole['M3']*g.ndata['T3'].reshape(n_frame,n_atom,3,3,3)).sum(-1).sum(-1).sum(-1).sum(-1))[:,None]/6
                    force_ext = force_ext - (g.ndata['T4'].reshape(-1,n_atom,3,3,3,3)*multipole['M3'][...,None]).sum(2).sum(2).sum(2)/6
            else:
                energy_ext = 0

            results = {}
            f = grad(
                energy + energy_ext,
                g.ndata['xyz'],
                grad_outputs=torch.ones_like(energy),
                create_graph=True,
                retain_graph=True,
            )[0] * CONV_UNIT_LEN  #* self.std_e

            if (len(T_keys)>0): 
                ft = grad(
                     energy + energy_ext,
                     [g.ndata[key] for key in T_keys],
                     grad_outputs=torch.ones_like(energy),
                     create_graph=True,
                     retain_graph=True)
                
                k = 0
                for key in T_keys:
                    results[key+'_gradient'] = ft[k]

                    if key == 'T0' and len(ft[k].shape)==1:
                        ft_qm = ft[k][:,None]*g.ndata['T1']
                    if key == 'T1' and len(ft[k].shape)==2:
                        ft_qm = ft_qm + (ft[k][...,None]*g.ndata['T2']).sum(1)
                    if key == 'T2' and len(ft[k].shape)==3:
                        ft_qm = ft_qm + (ft[k][...,None]*g.ndata['T3']).sum(1).sum(1)
                    if key == 'T3' and len(ft[k].shape)==4:
                        ft_qm = ft_qm + (ft[k][...,None]*g.ndata['T4']).sum(1).sum(1).sum(1)

                    k = k + 1
            
            results['qm_energy'] = energy #+ energy_ext
            results['ext_energy'] = energy_ext

            if self.reference_model:
                results['ext_force'] = force_ext.reshape(-1,n_atom,3)
            else:
                results['ext_force'] = 0

            results['qm_force'] = (-f + ft_qm).view(energy.shape[0],-1,3) - results['ext_force'] 
            results['qm_charge'] = charge
            results['feature_vector'] = feat
            results.update(multipole)

            if graph_rep is not None:
                results['graph_feature'] = graph_rep

        return results

    def training_step(self, batch, batch_id):
        g_qmmm, label, cell = batch

        if self.reference_model:
            if 'ref_force' not in label.keys() or 'ref_energy' not in label.keys():
                if self.train_total_energy:
                    Warning('qm_energy and qm_force should be the total one')
                else:
                    Warning('qm_energy and qm_force should be the difference between total one and reference model')
                # raise Exception('there should be ref_force when reference model is set to true')

        g_qmmm.ndata['xyz']['qm'].requires_grad_()
        if len(g_qmmm.ntypes)>1:
            g = dgl.node_type_subgraph(g_qmmm,['qm'])
            g.set_batch_num_edges({('qm','Hqm','qm'):g_qmmm.batch_num_edges(etype='Hqm'), 
                                   ('qm','Hqm2','qm'):g_qmmm.batch_num_edges(etype='Hqm2')})
            g.set_batch_num_nodes(g_qmmm.batch_num_nodes(ntype='qm'))

            if 'Tij1' not in g_qmmm.edges['Hqmmm'].data.keys():
                g_qmmm.ndata['xyz']['mm'].requires_grad_()
                T_tensor = compute_interacton_tensor(dgl.edge_type_subgraph(g_qmmm,[('mm','Hqmmm','qm')]))
                for key in T_tensor.keys():
                    g.ndata[key] = T_tensor[key] + g.ndata[key]
        else:
            g = g_qmmm

        bz = g.batch_size
        results = self(g, cell = cell)
        if 'Tij1' in g_qmmm.edges['Hqmmm'].data.keys():
            results['mm_force'] = compute_mm_analytical_grad(results,g_qmmm)
        else:
            fmm = grad(
                  results['qm_energy'] + results['ext_energy'],
                  g_qmmm.ndata['xyz']['mm'],
                  grad_outputs=torch.ones_like(results['qm_energy']),
                  create_graph=True,
                  retain_graph=True )[0]*CONV_UNIT_LEN
            results['mm_force'] = -fmm

        if self.out_model.name != "gp":
            if not self.train_total_energy and self.reference_model:
                loss_ene = torch.mean(abs(results['qm_energy'].squeeze() - label['qm_energy'].squeeze())**2)
            elif self.train_total_energy and 'ref_energy' in label.keys():
                loss_ene = torch.mean(abs((results['qm_energy']+results['ext_energy']).squeeze() - (label['qm_energy'] + label['ref_energy']).squeeze())**2)
            else:
                raise Exception('options are not compatible')
                #loss_ene = torch.mean(abs(results['qm_energy'].squeeze() - label['qm_energy'].squeeze())**2)
                #raise Exception('Error: option is not compatible')

            if len(g_qmmm.ntypes)>1:
                loss_mm_force = (abs(g_qmmm.ndata['f']['mm'] - results['mm_force'])**2).mean()
            else:
                loss_mm_force = 0

            if 'qm_force' in label.keys():
                if 'ref_force' in label.keys() and self.reference_model:
                    label_force = label['qm_force'] - label['ref_force']
                else:
                    label_force = label['qm_force'] 

                loss_force = torch.mean(abs(results['qm_force'] - label_force)**2) + loss_mm_force 
            else:
                if 'ref_force' in label.keys() and not self.train_total_energy:
                    label_force = g_qmmm.ndata['f']['qm'] - label['ref_force'].reshape(g_qmmm.ndata['f']['qm'].shape)
                else:
                    label_force = g_qmmm.ndata['f']['qm'] 

                if self.train_total_energy:
                    loss_force = torch.mean(abs((results['qm_force']+results['ext_force']).reshape(g_qmmm.ndata['f']['qm'].shape) - label_force)**2) + loss_mm_force 
                else:
                    loss_force = torch.mean(abs(results['qm_force'].reshape(g_qmmm.ndata['f']['qm'].shape) - label_force)**2) + loss_mm_force 

            if not self.mbis:
                if 'qm_charge' in label.keys():
                    loss_elec = torch.mean((abs(results['qm_charge'].view(label['qm_charge'].shape)-label['qm_charge']))**2)
                else:
                    loss_elec = torch.Tensor([0], device=loss_ene.device)
            else:
                loss_elec = 0
                if 'M0' in label.keys():
                    loss_elec = torch.mean((abs(results['qm_charge'].view(label['M0'].shape)-label['M0']))**2)
                if 'M1' in label.keys():
                    loss_elec = loss_elec + torch.mean((abs(results['M1'].view(label['M1'].shape)-label['M1']))**2)
                if 'M2' in label.keys():
                    loss_elec = loss_elec + 0.5*torch.mean((abs(results['M2'].view(label['M2'].shape)-label['M2'])*self.mask_M2[None,None])**2)
                if 'M3' in label.keys():
                    print('label_M3')
                    loss_elec = loss_elec + torch.mean((abs(results['M3'].view(label['M3'].shape)-label['M3'])*self.mask_M3[None,None])**2)/6

            ## L2 Loss
            l2_loss = Variable(torch.zeros(1,dtype=torch.float32,device=loss_elec.device), requires_grad=True)
            for name, param in self.named_parameters():
                if 'weight' in name:
                   l2_loss = l2_loss+(0.5 * torch.sum(torch.pow(param, 2)))

            L2_loss = self.wl2 * l2_loss
            loss = loss_ene + self.wf * loss_force + self.wq * loss_elec + L2_loss 

            self.log("averaged square deviation of energy", loss_ene*(CONV_UNIT_ENE**2))
            self.log("averaged square deviation of qm force", loss_force*(CONV_UNIT_FORCE**2))
            self.log("averaged square deviation of mm force", loss_mm_force*(CONV_UNIT_FORCE**2))
            self.log("averaged square deviation of multipole", loss_elec)
        return loss

    @torch.enable_grad()
    def validation_step(self, batch, batch_idx):

        # self.train()
        g_qmmm, label, cell = batch
        n_frame = g_qmmm.batch_size
        if len(g_qmmm.ntypes)>1:
            g = dgl.node_type_subgraph(g_qmmm,['qm'])
            g.set_batch_num_edges(g_qmmm.batch_num_edges(etype='Hqm'))
            g.set_batch_num_nodes(g_qmmm.batch_num_nodes(ntype='qm'))

            g_qmmm.ndata['xyz']['qm'].requires_grad_()
            if 'Tij1' not in g_qmmm.edges['Hqmmm'].data.keys():
                g_qmmm.ndata['xyz']['mm'].requires_grad_()
                T_tensor = compute_interacton_tensor(dgl.edge_type_subgraph(g_qmmm,[('mm','Hqmmm','qm')]),reference_model=self.reference_model)
                for key in T_tensor.keys():
                    g.ndata[key] = T_tensor[key] + g.ndata[key]
        else:
            g = g_qmmm

        n_atoms = label['qm_charge'].shape[1]
        results = self(g, cell = cell)

        if 'Tij1' in g_qmmm.edges['Hqmmm'].data.keys():
            results['mm_force'] = compute_mm_analytical_grad(results,g_qmmm)
        else:
            fmm = grad(
                  results['qm_energy']+ results['ext_energy'],
                  g_qmmm.ndata['xyz']['mm'],
                  grad_outputs=torch.ones_like(results['qm_energy']),
                  create_graph=True,
                  retain_graph=True )[0]*CONV_UNIT_LEN
            results['mm_force'] = -fmm

        if 'ref_energy' in label.keys() and self.reference_model:
            label['qm_energy'] = label['qm_energy'] - label['ref_energy']
            if 'qm_force' in label.keys():
                label['qm_force'] = label['qm_force'] - label['ref_force']

        metric = {}
        metric['amse_qmene']= abs(results['qm_energy'].squeeze() - label['qm_energy'].squeeze())*CONV_UNIT_ENE
        if 'qm_force' in label.keys():
            metric['amse_qmforce']= (abs(results['qm_force'] - label['qm_force']))*CONV_UNIT_FORCE
        else:
            if 'ref_force' in label.keys() and self.reference_model:
                label_force = g_qmmm.ndata['f']['qm'] - label['ref_force'].reshape(g_qmmm.ndata['f']['qm'].shape)
            else:
                label_force = g_qmmm.ndata['f']['qm'] 
            metric['amse_qmforce']= (abs(results['qm_force'] - label_force))*CONV_UNIT_FORCE

        if 'qm_charge' in label.keys():
            loss_charge = torch.mean((abs(results['qm_charge'].view(label['qm_charge'].shape)-label['qm_charge']))**2)

        if self.reference_model:
            if 'M0' in label.keys():
                energy_ext_ref = self.pooling1(g,g.ndata['T0']*label['M0'].view(-1))
                force_ext_ref = (g.ndata['T1'].reshape(-1,n_atoms,3)*label['M0'][...,None])
            else:
                energy_ext_ref = self.pooling1(g,g.ndata['T0']*label['qm_charge'].view(-1))
                force_ext_ref = (g.ndata['T1'].reshape(-1,n_atoms,3)*((label['qm_charge']).squeeze())[...,None])

            if 'M1' in label.keys() and 'T1' in g.ndata.keys():
                energy_ext_ref =  energy_ext_ref - ((label['M1']*g.ndata['T1'].reshape(n_frame,n_atoms,3)).sum(-1).sum(-1))
                force_ext_ref = force_ext_ref - (g.ndata['T2'].reshape(-1,n_atoms,3,3)*label['M1'][...,None]).sum(2)

            if 'M2' in label.keys() and 'T2' in g.ndata.keys():
                energy_ext_ref =  energy_ext_ref + 0.5*((label['M2']*g.ndata['T2'].reshape(n_frame,n_atoms,3,3)).sum(-1).sum(-1).sum(-1))
                force_ext_ref = force_ext_ref + (g.ndata['T3'].reshape(-1,n_atoms,3,3,3)*label['M2'][...,None]).sum(2).sum(2)/2

            if 'M3' in label.keys() and 'T3' in g.ndata.keys():
                energy_ext_ref =  energy_ext_ref - ((label['M3']*g.ndata['T3'].reshape(n_frame,n_atoms,3,3,3)).sum(-1).sum(-1).sum(-1).sum(-1))/6
                force_ext_ref = force_ext_ref - (g.ndata['T4'].reshape(-1,n_atoms,3,3,3,3)*label['M3'][...,None]).sum(2).sum(2).sum(2)/6

            metric['amse_qmmmene'] = abs(results['ext_energy']- energy_ext_ref)*CONV_UNIT_ENE
            metric['amse_qmmmforce'] = (abs(results['ext_force']- force_ext_ref).mean(-1)*CONV_UNIT_FORCE).mean(-1)

            qm_ene_total = results['ext_energy'].squeeze() + results['qm_energy'].squeeze()
            qm_ene_total_ref = label['qm_energy'] + energy_ext_ref
            qm_force_total = results['ext_force'] + results['qm_force']
            if 'qm_force' in label.keys():
                qm_force_total_ref = force_ext_ref + label['qm_force']
            else:
                qm_force_total_ref = force_ext_ref + g_qmmm.ndata['f']['qm'].reshape(-1,n_atoms,3)

            metric['amse_totene'] = abs(qm_ene_total - qm_ene_total_ref)*CONV_UNIT_ENE
            metric['amse_totforce'] = abs(qm_force_total - qm_force_total_ref).mean(-1).mean(-1)*CONV_UNIT_FORCE
        
            metric['qmmm_force_ref'] = force_ext_ref
            metric['qmmm_energy_ref'] = energy_ext_ref

        for key in results.keys():
            if results[key] is not None:
                results[key] = results[key].detach().cpu().numpy()

        for key in metric.keys():
            metric[key] = metric[key].detach().cpu().numpy()

        self.metric_key = list(metric.keys())
        self.results_key = list(results.keys())
        results['batch_idx'] = batch_idx
        return tuple(results.values()), tuple(metric.values()) 
    
    def validation_epoch_end(self, validation_step_outputs):
        pred, metric = map(list,zip(*validation_step_outputs))
        pred = list(zip(*pred))
        metric = list(zip(*metric))

        metrics = {}
        for i in range(len(self.metric_key)):
            metrics[self.metric_key[i]] = np.concatenate(metric[i])

        # preds = {}
        # for i in range(len(self.results_key)):
            # preds[self.results_key[i]] = np.concatenate(pred[i])

        if torch.distributed.is_initialized:
            self.ddp = True
            rank = torch.distributed.get_rank()
        else:
            rank = 0

        write_data(os.path.join(self.trainer.log_dir,f"metrics_{rank}.h5"),metrics)
        write_data(os.path.join(self.trainer.log_dir,f"pred_{rank}.h5"),pred)


    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        g, _ , cell = batch
        with torch.enable_grad():
            results = self(g, cell = cell)
        return results

    def configure_optimizers(self):
        return NotImplemented
