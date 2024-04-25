from distutils.log import error
import torch
import numpy as np
import dgl
from mlmm.baal.myheuristic import ForceVariance
from mlmm.utils.utils import check_value
import math
from mlmm.dataset.mydataset_alphachem import CONV_UNIT_ENE

class SamplerBase:
    def __init__(self, model=None, heuristic = ForceVariance):
        self.model = model
        self.heuristic = heuristic()

    def setup(self,model):
        self.model = model

    def run(self,ini_points):
        raise NotImplementedError


class OptimizerSampler(SamplerBase):
    def __init__(self, 
                 query_size = 10,
                 energy_alpha=10,
                 lr=0.001,
                 steps=200,
                 ref = 0,
                 sigma = 0.02,
                 scale = 0.1,
                 t1 = 0.0006,
                 t2 = 0.001,
                 wkernel = 1e-4,
                 increment = 10.0,
                 clip_value = 0.0,
                 shift = 1.8,
                 similarity_threhold= 0.05):
        super().__init__()
        self.lr = lr
        self.steps = steps
        self.energy_alpha = energy_alpha
        self.alpha = math.log(1e-20)/(-((energy_alpha/CONV_UNIT_ENE)**2))
        self.ref = ref
        self.sigma = sigma
        self.scale = scale
        self.t1 = t1
        self.t2 = t2
        self.wkernel = wkernel
        self.increment = increment
        self.clip_value = clip_value
        self.query_size = query_size
        self.shift = shift
        self.similarity_threhold = similarity_threhold

    def loss_fn(self, batch, results, feature = None):
        # qmmm_energy = results['qm_charge'] * batch[0].ndata['extfield_st']
        # energy = results['qm_energy'] + qmmm_energy
        bz = batch[0].batch_size
        natoms = batch[0].number_of_nodes()//bz
        qmmm_energy = (results[:,-natoms:] * batch[0].ndata['extfield_st'].reshape(bz,natoms)[...,None]).sum(1)
        energy = (results[:,0] + qmmm_energy)
        # energy = results[:,0] 
        prefactor = torch.exp(-self.alpha*(((energy-self.ref+self.shift)**2).mean(-1)))
        uncertainty = self.heuristic.get_uncertainties(results[:,1:1+natoms*3])
        #loss1 = prefactor*uncertainty
        # loss1 = -loss1.mean()
        loss1 = -uncertainty

        if self.wkernel > 0:
            if feature['graph_feature'] is not None:
                feat = feature['graph_feature']
                dis = (torch.norm(feat[:,None] - feat[None],dim=-1))
            else:
                feat = feature['feature_vector'].reshape(bz,natoms,-1)
                dis = (torch.norm(feat[:,None] - feat[None], dim = -1)).mean(-1)
            bandwidth = dis.detach().cpu().numpy()
            bandwidth = bandwidth[~np.eye(bandwidth.shape[0],dtype=bool)].reshape(bandwidth.shape[0],-1) 
            bandwidth = (torch.from_numpy(np.median(bandwidth,axis=-1)).to(feat.device))**2
            loss2 = torch.exp(-(dis**2)/bandwidth[:,None]).mean(-1)
        else:
            loss2 = 0
            dis = loss1

        loss = (prefactor*(loss1 + self.wkernel*loss2)).mean()
        return loss, uncertainty, [loss1.mean(),dis.mean()], prefactor, energy
    
    def update_alpha(self):
        self.energy_alpha = self.energy_alpha + self.increment
        self.alpha = math.log(1e-20)/(-((self.energy_alpha/CONV_UNIT_ENE)**2))
    
    def similarity_judge(self,feat,uncertainty):
        idx = torch.arange(feat.shape[0])
        idx3 = []
        while True:
            feat2 = feat[idx]
            dij = torch.norm(feat2[:,None] - feat2[None],dim=-1).mean(-1)
            dij = dij + torch.diag(torch.diag(dij)+0.2)
            idi, idj = torch.where(dij<self.similarity_threhold)
            if idi.shape[0] != 0:
                idxtmp = torch.cat([idx[idi],idx[idj]]).unique().cpu().numpy()
                idx3.append(np.setdiff1d(idx,idxtmp))
                idx = torch.where((uncertainty[idx[idi]]>uncertainty[idx[idj]]).cpu(),idx[idi],idx[idj])
                idx = idx.unique()
            else:
                break
        
        if idx.shape[0] != feat.shape[0]: 
            idx3.append(idx)
            idx = np.concatenate(idx3)

        return idx

        
    def run(self, batch):
        # for var in batch:
            # if isinstance(var,dgl.heterograph):
                # has_input = True
        
        # if not has_input:
            # raise TypeError('no dgl input')
        bz = batch[0].batch_size
        natoms = batch[0].number_of_nodes()//bz
        if torch.rand([1]) > 0.5:
            scale_factor = 1 + (torch.rand([1]) - 0.5)*2*self.scale
            batch[0].ndata['extfield'] = batch[0].ndata['extfield']*scale_factor.to(batch[0].device)
            batch[0].ndata['extfield_st'] = batch[0].ndata['extfield_st']*scale_factor.to(batch[0].device)
 
            std1 = batch[0].ndata['extfield'].std()
            std2 = batch[0].ndata['extfield_st'].std()
            batch[0].ndata['xyz'] = batch[0].ndata['xyz'] + self.sigma*torch.randn(batch[0].ndata['xyz'].shape).to(batch[0].device)
            batch[0].ndata['extfield'] = batch[0].ndata['extfield'] + std1*self.sigma*torch.randn(batch[0].ndata['extfield'].shape).to(batch[0].device)
            batch[0].ndata['extfield_st'] = batch[0].ndata['extfield_st'] + std2*self.sigma*torch.randn(batch[0].ndata['extfield_st'].shape).to(batch[0].device)
  
        batch[0].ndata['xyz'].requires_grad_()
        batch[0].ndata['extfield'].requires_grad_()
        batch[0].ndata['extfield_st'].requires_grad_()

        opt = torch.optim.Adam([batch[0].ndata['xyz'],batch[0].ndata['extfield'],batch[0].ndata['extfield_st']], lr=self.lr)

        i = 0
        while i < self.steps: 
            results, feature = self.model.predict_step(batch,0)#, batch_idx, dataloader_idx=dataloader_idx)
            # results= self.model.predict_step(batch,0)#, batch_idx, dataloader_idx=dataloader_idx)
            # results = torch.stack(results)
            loss, uncertainty, metric, prefactor,energy = self.loss_fn(batch,results,feature=feature)

            if i == 0:
                ini_loss = loss.detach().clone()
                ini_loss1 = metric[0].detach().clone()
                ini_loss2 = metric[1].detach().clone()
            #print(f"ini_loss:{ini_loss}")
            #print(f"ini_loss1:{ini_loss1}")
            #print(f"ini_loss2:{ini_loss2}")
            # loss = self.loss_fn(batch,results)

            opt.zero_grad()
            loss.backward()
            if self.clip_value > 0:
                torch.nn.utils.clip_grad_norm_([batch[0].ndata['xyz'],batch[0].ndata['extfield'],batch[0].ndata['extfield_st']],max_norm=self.clip_value)
            opt.step()
            print(f"step:{i} loss:{loss}, NSDE:{metric[0]}, kernel loss:{metric[1]}") 
            print(f"step:{i} alpha:{self.alpha}, uncertainty:{uncertainty}")
            print(f"energy:{energy.mean(-1)}")
            print(f"prefactor:{prefactor}")

            i = i + 1

        results, feature = self.model.predict_step(batch,0)
        loss, uncertainty, metric, prefactor, energy = self.loss_fn(batch,results,feature=feature)
        idx = torch.where(uncertainty>self.t1)[0]
       # idx2 = torch.where(results[:,0].mean(-1) <= self.energy_alpha/CONV_UNIT_ENE + self.ref - self.shift)[0]
        #uncertainty_std = uncertainty.std()
        #idx2 = torch.where(uncertainty - uncertainty.mean() < 3.5*uncertainty_std)[0]
        uncertainty_median = torch.median(uncertainty)
        prefactor_median = torch.median(prefactor)
        idx2 = torch.where( uncertainty/uncertainty_median < 3.0)[0]
        #idx3 = torch.where( prefactor/prefactor_median < 10000)[0]
        #idx2 = np.intersect1d(idx2.cpu().numpy(), idx3.cpu().numpy())
        idx = np.intersect1d(idx.cpu().numpy(),idx2.cpu().numpy())
        #idx = np.intersect1d(idx.cpu().numpy(),idx2)
        print(f"step:{i} loss:{loss}, NSDE:{metric[0]}, kernel loss:{metric[1]}") 
        print(f"step:{i} alpha:{self.alpha}, uncertainty:{uncertainty}")
        print(f"energy:{energy.mean(-1)}")
        print(f"prefactor:{prefactor}")

        if uncertainty.mean() < self.t2 and self.energy_alpha < 1500:
            self.update_alpha()
            # self.energy_alpha = self.energy_alpha + self.increment
            # self.alpha = math.log(0.5)/(-((self.energy_alpha/CONV_UNIT_ENE)**2))
            # self.alpha = self.alpha - self.increment
            # if self.alpha < 0:
                # self.alpha = 0
        
        if idx.shape[0] > self.query_size:
            feat = feature['feature_vector'].reshape(bz,natoms,-1)
            idx3 = self.similarity_judge(feat[idx],uncertainty[idx])
            idx = idx[idx3]
            _, idx3 = torch.sort(uncertainty[idx],descending=True)
            idx = idx[idx3.cpu()]
            idx = idx[:self.query_size]
        print(f"meet condition: {idx.shape[0]}")
        print(f"select uncertainty: {uncertainty[idx]}")

        g = []
        for i in range(idx.shape[0]):
            g.append(dgl.slice_batch(batch[0], int(idx[i])))
        g = dgl.batch(g)

        out = {}
        out['qm_energy'] = results[idx,0].mean(-1).detach().cpu()
        out['qm_force'] = (results[idx,1:natoms*3+1].mean(-1)).reshape(idx.shape[0],natoms,3).detach().cpu()
        out['qm_charge'] = results[idx,natoms*3+1:].mean(-1).detach().cpu()
        return g, out
        

