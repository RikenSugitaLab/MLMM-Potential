import torch
from torch import nn
from torch.nn import init

from mlmm.nn.activations import shifted_softplus, swish
from mlmm.nn.base import Dense
from mlmm.nn.initializers import zeros_initializer
from typing import Dict, Union, List


__all__ = ["MLP", "TiledMultiLayerNN", "ElementalGate"]


class MLP(nn.Module):
    """Multiple layer fully connected perceptron neural network.

    Args:
        n_in (int): number of input nodes.
        n_out (int): number of output nodes.
        n_hidden (list of int or int, optional): number hidden layer nodes.
            If an integer, same number of node is used for all hidden layers resulting
            in a rectangular network.
            If None, the number of neurons is divided by two after each layer starting
            n_in resulting in a pyramidal network.
        n_layers (int, optional): number of layers.
        activation (ActBaseClass, optional): activation function. All hidden layers would
            the same activation function except the output layer that does not apply
            any activation function.
        dropout (bool): if adding dropout layer
        p (float): dropping rate

    """

    def __init__(
        self, n_in: int, n_out: int, n_hidden=None, n_layers=2, activation = None, 
        dropout = False, p = 0.1
    ):
        super(MLP, self).__init__()
        # get list of number of nodes in input, hidden & output layers
        if n_hidden is None:
            c_neurons = n_in
            self.n_neurons = []
            for i in range(n_layers):
                self.n_neurons.append(c_neurons)
                c_neurons = c_neurons // 2
            self.n_neurons.append(n_out)
        else:
            # get list of number of nodes hidden layers
            if type(n_hidden) is int:
                n_hidden = [n_hidden] * (n_layers - 1)
            self.n_neurons = [n_in] + n_hidden + [n_out]

        # assign a Dense layer (with activation function) to each hidden layer
        layers = []
        for i in range(n_layers - 1):
            layers.append( Dense(self.n_neurons[i], self.n_neurons[i + 1], activation=activation))
            if dropout:
                layers.append(torch.nn.Dropout(p=p))
        # layers = [
            # Dense(self.n_neurons[i], self.n_neurons[i + 1], activation=activation)
            # for i in range(n_layers - 1)
        # ]
        # assign a Dense layer (without activation function) to the output layer
        layers.append(Dense(self.n_neurons[-2], self.n_neurons[-1], activation=None))
        # put all layers together to make the network
        self.out_net = nn.Sequential(*layers)

    def forward(self, inputs):
        """Compute neural network output.

        Args:
            inputs (torch.Tensor): network input.

        Returns:
            torch.Tensor: network output.

        """
        return self.out_net(inputs)


class TiledMultiLayerNN(nn.Module):
    """
    Tiled multilayer networks which are applied to the input and produce n_tiled different outputs.
    These outputs are then stacked and returned. Used e.g. to construct element-dependent prediction
    networks of the Behler-Parrinello type.

    Args:
        n_in (int): number of input nodes
        n_out (int): number of output nodes
        n_tiles (int): number of networks to be tiled
        n_hidden (int): number of nodes in hidden nn (default 50)
        n_layers (int): number of layers (default: 3)
    """

    def __init__(
        self, n_in, n_out, n_tiles, n_hidden=50, n_layers=3, activation=shifted_softplus
    ):
        super(TiledMultiLayerNN, self).__init__()
        self.mlps = nn.ModuleList(
            [
                MLP(
                    n_in,
                    n_out,
                    n_hidden=n_hidden,
                    n_layers=n_layers,
                    activation=activation,
                )
                for _ in range(n_tiles)
            ]
        )

    def forward(self, inputs):
        """
        Args:
            inputs (torch.Tensor): Network inputs.

        Returns:
            torch.Tensor: Tiled network outputs.

        """
        return torch.cat([net(inputs) for net in self.mlps], 2)


class ElementalGate(nn.Module):
    """
    Produces a Nbatch x Natoms x Nelem mask depending on the nuclear charges passed as an argument.
    If onehot is set, mask is one-hot mask, else a random embedding is used.
    If the trainable flag is set to true, the gate values can be adapted during training.

    Args:
        elements (set of int): Set of atomic number present in the data
        onehot (bool): Use one hit encoding for elemental gate. If set to False, random embedding is used instead.
        trainable (bool): If set to true, gate can be learned during training (default False)
    """

    def __init__(self, elements, onehot=True, trainable=False):
        super(ElementalGate, self).__init__()
        self.trainable = trainable

        # Get the number of elements, as well as the highest nuclear charge to use in the embedding vector
        self.nelems = len(elements)
        maxelem = int(max(elements) + 1)

        self.gate = nn.Embedding(maxelem, self.nelems)

        # if requested, initialize as one hot gate for all elements
        if onehot:
            weights = torch.zeros(maxelem, self.nelems)
            for idx, Z in enumerate(elements):
                weights[Z, idx] = 1.0
            self.gate.weight.data = weights

        # Set trainable flag
        if not trainable:
            self.gate.weight.requires_grad = False

    def forward(self, atomic_numbers):
        """
        Args:
            atomic_numbers (torch.Tensor): Tensor containing atomic numbers of each atom.

        Returns:
            torch.Tensor: One-hot vector which is one at the position of the element and zero otherwise.

        """
        return self.gate(atomic_numbers)
