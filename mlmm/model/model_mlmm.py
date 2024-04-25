import torch
import dgl
from torch import nn
from mlmm.nn.base import Dense
from mlmm.nn.blocks import MLP
from torch.autograd import grad
from dgl.nn.pytorch.glob import SumPooling

__all__=["model_mlmm","general_derivative"]

def general_derivative(fx, dx, hessian=False, create_graph=True, retain_graph=True):
    fx_shape = fx.size()
    dx_shape = dx.size()

    # print('=======================')
    # print(fx_shape, 'FX')
    # print(dx_shape, 'DX')

    fx = fx.view(fx_shape[0], -1)
    # print(fx.shape, 'REFORM')

    dfxdx = torch.stack(
        [grad(fx[..., i], dx, torch.ones_like(fx[:, i]), create_graph=create_graph, retain_graph=retain_graph)[0] for
         i in range(fx.size(1))], dim=1
    )

    if not hessian:
        dfxdx = dfxdx.view(*fx_shape, *dx_shape[1:])
        # TODO: squeeze for energies might be dangerous
        dfxdx = dfxdx.squeeze(dim=1)  # Saver squeeze for batch size of 1
    else:
        dfxdx = dfxdx.view(fx_shape[0], fx_shape[-1] * fx_shape[-2], -1)

    # print(dfxdx.shape, 'DFXDX')
    # print('=======================')

    return dfxdx

class model_mlmm(nn.Module):
    def __init__(self, 
                 rep_model, 
                 n_hidden=None,
                 n_layers=2,
                 mean_e=0,
                 std_e=1):

        super(model_mlmm,self).__init__()

        self.rep_model=rep_model
        self.pooling = SumPooling()
        self.mean_e  = mean_e
        self.std_e = std_e

        self.oute = MLP(self.rep_model.n_atom_basis,1, n_hidden=n_hidden, n_layers=n_layers)
        self.outq = MLP(self.rep_model.n_atom_basis,1, n_hidden=n_hidden, n_layers=n_layers)

    def forward(self, g, extfield, cell=None, total_charge=-1.0):

        g.ndata['xyz'].requires_grad_()
        feat, mu = self.rep_model(g, extfield, cell=cell)
        energy = self.pooling(g,self.oute(feat))
        n_atom = g.batch_num_nodes()[0]

    #    dipole = -general_derivative(E, extfield)
    #    dipole = torch.sum(dipole,1)
    #    d2EdRdF = general_derivative(-dipole, inputs[Properties.R])
    #    d2EdRdF = d2EdRdF.permute(0, 2, 3, 1)
        charge = self.outq(feat)
        charge = charge.reshape(-1,n_atom)
        charge = charge + (total_charge - charge.sum(-1))[:,None]/n_atom

        f = grad(
            energy,
            g.ndata['xyz'],
            grad_outputs=torch.ones_like(energy),
            create_graph=True,
            retain_graph=True,
        )[0] * self.std_e

        energy = (energy * self.std_e) + self.mean_e

        results = {}
        results['qm_energy'] = energy
        results['qm_force'] = -f.view(energy.shape[0],-1,3)
        results['qm_charge'] = charge
        print(results['qm_energy'].dtype)
        return results

