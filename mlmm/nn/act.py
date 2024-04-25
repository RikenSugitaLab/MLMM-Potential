import torch 
from torch import nn
from torch.cuda.amp import autocast
from mlmm.nn import Dense
import numpy as np
import math
import os
if "use_fp16" in os.environ:
    if os.environ['use_fp16'] == 'True':
        use_fp16=True
    elif os.environ['use_fp16'] == 'False':
        use_fp16=False
    else:
        raise AssertionError("wrong setting, use_fp16 can only be True or False")
else:
    use_fp16=False

__all__ = ["AdaptiveComputationTime","_gen_timing_signal"]
class AdaptiveComputationTime(nn.Module):
    r"""AdaptiveComputationTime used in 3DT module.

    Args:
        n_in (int): number of input (i.e. atomic embedding) dimensions.
        
    """

    def __init__(
        self,
        n_atom_basis,
        activation=None,
        bias = 1.,
    ):
        super(AdaptiveComputationTime, self).__init__()
        
        self.threshold = 1. - 0.1
        self.init_bias= bias
        

#         ### v3:
#         self.halt_unit = Dense(2*n_atom_basis, 1, activation=None )        
#         nn.init.constant_(self.halt_unit.bias, self.init_bias)
        
        ### v4:
        self.halt_unit = Dense(3*n_atom_basis, 1, activation=None )        
        nn.init.constant_(self.halt_unit.bias, self.init_bias)
        
        self.sigma = nn.Sigmoid()

    @autocast(enabled=use_fp16)       
    def forward(self, g, state, e, fn, time_enc, max_hop, cell=None
                ):
        ### state:(N_b, N_a, atom_basis)   
        ### time_enc: [1, LEN, atom_basis)
        self.g = g.local_var()

        halting_probability = torch.zeros(state.size()[0],state.size()[1]).cuda() ### (N_b, N_a)  
        remainders = torch.zeros(state.size()[0],state.size()[1]).cuda() ### (N_b, N_a)  
        n_updates = torch.zeros(state.size()[0],state.size()[1]).cuda() ### (N_b, N_a)  
        previous_state = torch.zeros_like(state).cuda() ### (N_b, N_a, atom_basis)  
        step = 0

        while( ((halting_probability<self.threshold) & (n_updates<max_hop)).byte().any() ):
            t = time_enc[:, step, :].unsqueeze(1).type_as(state) ### [1,1,dim]
            t_batch = time_enc[:, step, :].repeat([e.shape[0],1])
            state_embed = torch.cat( [state, t_batch, e], axis=-1) ### Time embedding

            with autocast(enabled=False):
                p = self.sigma(self.halt_unit(state_embed.float()))### (N_b, N_a)
            

            ### Mask of inputs which have not halted yet:
            still_running = (halting_probability < 1.0).float()

            ### Mask of inputs which halted at this step:
            new_halted = (halting_probability + p*still_running > self.threshold).float() * still_running

            ### Mask of inputs which haven't halted, and did not hallt at this step:
            still_running = (halting_probability + p*still_running <= self.threshold).float()

            ### Add the halting probability for this step 
            ### to the halting probabilities for inputs which haven't halted yet:
            halting_probability = halting_probability +  p*still_running

            ### Compute remainders for the inputs which halted at this step
            remainders = remainders + new_halted * (1. - halting_probability)

            ### Add the remainders to those which halted at this step:
            halting_probability = halting_probability + new_halted*remainders

            ### Increment n_updates for all inputs which are still running:
            n_updates = n_updates + still_running + new_halted

            ### Compute the weight to be applied to the new state and output
            ### 0 when halted
            ### p when not halted yet
            ### remainders when halted at this step:
            update_weights = p*still_running + new_halted*remainders            
            
            ###state = fn(e, state, t_fn, r_iij, neighbors, neighbor_mask_iij, f_ij=f_iij)
            
            ### Apply Transition on the state:
            t_fn = t.squeeze(1) ### [1,dim]
            state = fn(self.g, state, t_fn, cell=cell)

            ### update running part in the weighted states and keep the rest:
            previous_state = ((state * update_weights) + \
                              (previous_state * (1. - update_weights)))
            ### Notice that indeed we return the previous_state
            
#             ### Apply Transition on the state:
#             t_fn = t.squeeze(1) ### [1,dim]
#             if ((halting_probability<self.threshold) & (n_updates<max_hop)).byte().any():
#                 state = fn(e, state, t_fn, r_iij, neighbors, neighbor_mask_iij, f_ij=f_iij)

            step += 1

        return previous_state, (remainders,n_updates)               
                
@autocast(enabled=use_fp16)
def _gen_timing_signal(length, channels, min_timescale=1.0, max_timescale=1.0e4):
    """
    Generate a [1, lenght, channels] time embedding
    """
    
    position = np.arange(length)
    num_timescales = channels // 2
    log_timescale_increment = ( math.log(float(max_timescale) / float(min_timescale)) / (float(num_timescales)-1))
    inv_timescales = min_timescale * np.exp(np.arange(num_timescales).astype(np.float) * -log_timescale_increment)
    scaled_time = np.expand_dims(position,1) * np.expand_dims(inv_timescales,0)
    
    signal = np.concatenate([np.sin(scaled_time), np.cos(scaled_time)], axis=1)
    signal = np.pad(signal, [[0,0], [0,channels%2]],
                    'constant', constant_values=[0.0,0.0])
    signal = signal.reshape([1,length,channels])
    
    if use_fp16:
        return torch.from_numpy(signal).half()
    else:
        return torch.from_numpy(signal).type(torch.FloatTensor)

 
