import random
import numpy as np
import argparse
import torch
import dgl
import h5py
import os
from mlmm.nn import Dense

FALSY_STRINGS = {"off", "false", "0"}
TRUTHY_STRINGS = {"on", "true", "1"}
CONV_UNIT_LEN    = 0.52917721092
key_list = ['T0','T1','T2','T3']

# def set_seed(args=None):
    # seed = 1 if not args else args.seed
    # 
    # random.seed(seed)
    # np.random.seed(seed)
# 
    # torch.manual_seed(seed)
    # torch.cuda.manual_seed_all(seed)
    # torch.backends.cudnn.deterministic = True
    # torch.backends.cudnn.benchmark = False
# 
    # dgl.random.seed(seed)
# 
# 
def make_log_dir(args):
    # make and return
    if args.log_path=="":
        output_dir = os.path.join(args.model_path,f"log&check")
    else: 
        output_dir = os.path.join(args.model_path,args.log_path)

    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(output_dir,f"run_log")
    os.makedirs(log_dir, exist_ok=True)
    return log_dir,output_dir

def moment_update(model,model_ema,m):
    """ model_ema = m* model_ema + (1-m)model """
    for p1, p2 in zip(model.parameters(), model_ema.parameters()):
        p2.data.mul_(m).add_(1-m,p1.detach().data)

def dist_collect(x,eval=False):
    x = x.contiguous()
    out_list = [torch.zeros_like(x,device=x.device, dtype=x.dtype) for _ in range(torch.distributed.get_world_size())]
    torch.distributed.all_gather(out_list,x)
    if eval:
        world_size = torch.distributed.get_world_size()
        back_id = []
        for i in range(len(out_list)):
            try:
                if len(out_list[i].shape)==0:
                   size1 = 1
                   out_list[i] = out_list[i][None]
                else:
                   size1 = out_list[i].shape[0]
            except:
                print(out_list[i])
                continue
            back_id.append(torch.arange(size1)*world_size+i)
        back_id = torch.cat(back_id,dim=0).cuda()
        out_list = torch.cat(out_list,dim=0)
        out = torch.zeros_like(out_list) 
        out.index_copy_(0,back_id,out_list)
        return out
    else:
        return torch.cat(out_list,dim=0)

def write_data(filename,data):
    with h5py.File(filename,'w') as f:
        for i in data: 
           f[i]=data[i]

def read_filelist(filename):
    with open(filename,'r') as f:
        pos_file=f.readlines()
  
        for i in range(len(pos_file)):
           pos_file[i]=pos_file[i].strip('\n')
    return pos_file

def check_value(x):
    num = torch.isnan(x).sum() + torch.isinf(x).sum()
    if num>0:
        out = True
    else:
        out = False
    return out

def dict2namespace(config):
    namespace = argparse.Namespace()
    for key, value in config.items():
        if isinstance(value,dict):
            new_value = dict2namespace(value)
        else:
            new_value = value
        setattr(namespace,key,new_value)
    return namespace

def bool_flag(s):
    """
    Parse boolean arguments from the command line.
    """
    if s.lower() in FALSY_STRINGS:
        return False
    elif s.lower() in TRUTHY_STRINGS:
        return True
    else:
        raise argparse.ArgumentTypeError("invalid value for a boolean flag")

def add_sn(m):
    for name, layer in m.named_children():
         m.add_module(name, add_sn(layer))
    if isinstance(m, (Dense)):
         return torch.nn.utils.spectral_norm(m)
    else:
         return m
    
def compute_interacton_tensor(g1, q=None, gamma=0, delta = 0):
    g = g1.local_var()
    if len(g.ntypes)>1:
        dis_vec = dgl.ops.v_sub_u(g,g.ndata['xyz']['qm'].detach(),g.ndata['xyz']['mm'])/CONV_UNIT_LEN
    else:
        if q is not None:
            g.ndata['c'] = q
        dis_vec = dgl.ops.v_sub_u(g,g.ndata['xyz'],g.ndata['xyz'])/CONV_UNIT_LEN

    dis = torch.norm(dis_vec+1e-8, dim=-1)
    results = {}
    if len(g.ntypes)>1:
        results['T0'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,1/dis,g.ndata['c']['mm']))
        results['T1'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,dis_vec*torch.pow(dis,-3)[...,None],g.ndata['c']['mm'][...,None]))
    else:
        coef = 0.5*torch.erf((dis - delta)*gamma) + 0.5
        results['T0'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,coef/dis,g.ndata['c']))
        results['T1'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,dis_vec*(coef*torch.pow(dis,-3))[...,None],g.ndata['c'][...,None]))

    I = torch.eye(3,device=dis.device)
    if 'T2' in g.ndata.keys():
        if len(g.ntypes)>1:
            results['T2'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,(3*dis_vec[...,None]*dis_vec[:,None] - (dis**2)[:,None,None]*I[None])/torch.pow(dis[:,None,None],5),g.ndata['c']['mm'][...,None,None]))
        else:
            results['T2'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,coef[:,None,None]*(3*dis_vec[...,None]*dis_vec[:,None] - (dis**2)[:,None,None]*I[None])/torch.pow(dis[:,None,None],5),g.ndata['c'][...,None,None]))

    if 'T3' in g.ndata.keys():
        results['T3'] = 15*dis_vec[:,:,None,None]*dis_vec[:,None,:,None]*dis_vec[:,None,None,:] - \
                        3*(dis**2)[:,None,None,None]*(dis_vec[:,:,None,None]*I[None,None,:,:] + \
                                                      dis_vec[:,None,:,None]*I[None,:,None,:] + \
                                                      dis_vec[:,None,None,:]*I[None,:,:,None])
        
        if len(g.ntypes)>1:
            results['T3'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,results['T3']/torch.pow(dis[:,None,None,None],7),g.ndata['c']['mm'][...,None,None,None]))
        else:
            results['T3'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,coef[:,None,None,None]*results['T3']/torch.pow(dis[:,None,None,None],7),g.ndata['c'][...,None,None,None]))

    if 'T4' in g.ndata.keys():
        results['T4'] = 105*dis_vec[:,:,None,None,None]*dis_vec[:,None,:,None,None]*dis_vec[:,None,None,:,None]*dis_vec[:,None,None,None,:] - \
                        15*(dis**2)[:,None,None,None,None]*(dis_vec[:,:,None,None,None]*dis_vec[:,None,:,None,None]*I[None,None,None,:,:] + \
                                                            dis_vec[:,:,None,None,None]*dis_vec[:,None,None,:,None]*I[None,None,:,None,:] + \
                                                            dis_vec[:,:,None,None,None]*dis_vec[:,None,None,None,:]*I[None,None,:,:,None] + \
                                                            dis_vec[:,None,:,None,None]*dis_vec[:,None,None,:,None]*I[None,:,None,None,:] + \
                                                            dis_vec[:,None,:,None,None]*dis_vec[:,None,None,None,:]*I[None,:,None,:,None] + \
                                                            dis_vec[:,None,None,:,None]*dis_vec[:,None,None,None,:]*I[None,:,:,None,None]) 

        results['T4'] = results['T4'] + 3*torch.pow(dis,4)[:,None,None,None,None]*(I[None,:,:,None,None]*I[None,None,None,:,:] + I[None,:,None,:,None]*I[None,None,:,None,:] + I[None,:,None,None,:]*I[None,None,:,:,None]) 
        if len(g.ntypes)>1:
            results['T4'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,results['T4']/torch.pow(dis[:,None,None,None,None],9),g.ndata['c']['mm'][...,None,None,None,None]))
        else:
            results['T4'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,coef[:,None,None,None,None]*results['T4']/torch.pow(dis[:,None,None,None,None],9),g.ndata['c'][...,None,None,None,None]))
    return results

def compute_mm_analytical_grad(results,g):
    g_qmmm = g.edge_type_subgraph([('mm','Hqmmm','qm')])
    grad = -dgl.ops.copy_e_sum(g_qmmm.reverse(),dgl.ops.v_mul_e(g_qmmm,results['T0_gradient'][:,None],g.edata['Tij1'][('mm', 'Hqmmm', 'qm')]))
    grad = grad - dgl.ops.copy_e_sum(g_qmmm.reverse(),dgl.ops.e_mul_v(g_qmmm,g.edata['Tij2'][('mm', 'Hqmmm', 'qm')],results['T1_gradient'][:,None]).sum(-1))

    if 'Tij3' in g.edata.keys():
        grad = grad - dgl.ops.copy_e_sum(g_qmmm.reverse(),dgl.ops.e_mul_v(g_qmmm,g.edata['Tij3'][('mm', 'Hqmmm', 'qm')],results['T2_gradient'][...,None]).sum(1).sum(1))
    if 'Tij4' in g.edata.keys():
        grad = grad - dgl.ops.copy_e_sum(g_qmmm.reverse(),dgl.ops.e_mul_v(g_qmmm,g.edata['Tij4'][('mm', 'Hqmmm', 'qm')],results['T3_gradient'][...,None]).sum(1).sum(1).sum(1))
    return grad

def compute_interacton_tensor2(g1, dis, dis_vec, q=None, gamma=0, delta=0):
    g = g1.local_var()

    if q is not None:
        g.ndata['c'] = q

    dis_vec = dis_vec/CONV_UNIT_LEN
    dis = dis/CONV_UNIT_LEN

    results = {}
    if len(g.ntypes)>1:
        results['T0'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,1/dis,g.ndata['c']['mm']))
        results['T1'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,dis_vec*torch.pow(dis,-3)[...,None],g.ndata['c']['mm'][...,None]))
    else:
        coef = 0.5*torch.erf((dis - delta)*gamma) + 0.5
        results['T0'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,coef/dis,g.ndata['c']))
        results['T1'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,dis_vec*(coef*torch.pow(dis,-3))[...,None],g.ndata['c'][...,None]))

    I = torch.eye(3,device=dis.device)
    if 'T2' in g.ndata.keys() and 'T3' in g.ndata.keys():
        if len(g.ntypes)>1:
            results['T2'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,(3*dis_vec[...,None]*dis_vec[:,None] - (dis**2)[:,None,None]*I[None])/torch.pow(dis[:,None,None],5),g.ndata['c']['mm'][...,None,None]))
        else:
            results['T2'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,coef[:,None,None]*(3*dis_vec[...,None]*dis_vec[:,None] - (dis**2)[:,None,None]*I[None])/torch.pow(dis[:,None,None],5),g.ndata['c'][...,None,None]))

    if 'T3' in g.ndata.keys() and 'T4' in g.ndata.keys():
        results['T3'] = 15*dis_vec[:,:,None,None]*dis_vec[:,None,:,None]*dis_vec[:,None,None,:] - \
                        3*(dis**2)[:,None,None,None]*(dis_vec[:,:,None,None]*I[None,None,:,:] + \
                                                      dis_vec[:,None,:,None]*I[None,:,None,:] + \
                                                      dis_vec[:,None,None,:]*I[None,:,:,None])
        
        if len(g.ntypes)>1:
            results['T3'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,results['T3']/torch.pow(dis[:,None,None,None],7),g.ndata['c']['mm'][...,None,None,None]))
        else:
            results['T3'] = dgl.ops.copy_e_sum(g,dgl.ops.e_mul_u(g,coef[:,None,None,None]*results['T3']/torch.pow(dis[:,None,None,None],7),g.ndata['c'][...,None,None,None]))

    return results


