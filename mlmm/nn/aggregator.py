import numpy as np
import torch
import dgl
from mlmm.nn import Dense
__all__ = ['aggregator']
class aggregator(torch.nn.Module):
    def __init__(self,
                 mode: str = 'all') -> None:
        super().__init__()
        self.mode = mode
        if self.mode == 'all':
            self.combine = Dense(4,1)

    def forward(self,g,m): 
        if self.mode == 'sum' or self.mode == 'all':
            m_sum = dgl.ops.copy_e_sum(g,m)
            m2 = m_sum

        if self.mode == 'min' or self.mode == 'all':
            m_min = dgl.ops.copy_e_min(g,m)
            m2 = m_min

        if self.mode == 'max' or self.mode == 'all':
            m_max = dgl.ops.copy_e_max(g,m)
            m2 = m_max

        if self.mode == 'mean' or self.mode == 'all':
            m_mean = dgl.ops.copy_e_mean(g,m)
            m2 = m_mean

        if self.mode == 'all':
            m = torch.stack([m_sum,m_min,m_max,m_mean],dim=-1)
            m = torch.squeeze(self.combine(m))
            return m
        else:
            return m2

