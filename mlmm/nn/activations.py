import numpy as np
from torch.nn import functional
import torch

__all__ = ["shifted_softplus","swish","swish2","SiLU","Swish"]

def shifted_softplus(x):
    r"""Compute shifted soft-plus activation function.

       .. math::
          y = \ln\left(1 + e^{-x}\right) - \ln(2)

       Args:
           x (torch.Tensor): input tensor.

       Returns:
           torch.Tensor: shifted soft-plus of input.

       """
    return functional.softplus(x) - np.log(2.0)

def swish(x):
    #return x*functional.sigmoid(0.5*x)/1.1
    return x*functional.sigmoid(x)

def swish2(x):
    return x*functional.sigmoid(0.5*x)/1.1


class Swish(torch.nn.Module):
    def __init__(self,alpha,beta):
        super(Swish,self).__init__()
        self.alpha = alpha
        self.beta = beta
    
    def forward(self,x):
        y = x*functional.sigmoid(self.alpha*x)/self.beta
        return y

class SiLU(torch.nn.Module):
    def __init__(self,alpha,beta):
        super(SiLU,self).__init__()
        self.alpha = alpha
        self.beta = beta
    
    def forward(self,x):
        y = self.alpha*x/(1+torch.exp(-self.beta*x))
        return y
