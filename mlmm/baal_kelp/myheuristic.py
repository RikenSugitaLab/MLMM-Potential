import torch
from baal.active.heuristics import AbstractHeuristic
import types
import warnings
from collections.abc import Sequence
from functools import wraps as _wraps
from typing import List

import numpy as np
import scipy.stats
import torch
from scipy.special import xlogy
from torch import Tensor

from baal.utils.array_utils import to_prob
__all__=["MyAbstractHeuristic","ForceVariance"]

available_reductions = {
    "max": lambda x: torch.max(x,dim=-1).values,# axis=tuple(range(1, x.ndim))),
    "min": lambda x: torch.min(x, axis=tuple(range(1, x.ndim))),
    "mean": lambda x: torch.mean(x, axis=tuple(range(1, x.ndim))),
    "sum": lambda x: torch.sum(x, axis=tuple(range(1, x.ndim))),
    "none": lambda x: x,
}

def _shuffle_subset(data: torch.Tensor, shuffle_prop: float) -> torch.Tensor:
    to_shuffle = np.nonzero(np.random.rand(data.shape[0]) < shuffle_prop)[0]
    data[to_shuffle, ...] = data[np.random.permutation(to_shuffle), ...]

class MyAbstractHeuristic:
    """
    Abstract class that defines a Heuristic.

    Args:
        shuffle_prop (float): shuffle proportion.
        reverse (bool): True if the most uncertain sample has the highest value.
        reduction (Union[str, Callable]): Reduction used after computing the score.
    """

    def __init__(self, shuffle_prop=0.0, reverse=False, reduction="none"):
        self.shuffle_prop = shuffle_prop
        self.reversed = reverse
        assert reduction in available_reductions or callable(reduction)
        self._reduction_name = reduction
        self.reduction = reduction if callable(reduction) else available_reductions[reduction]

    def compute_score(self, predictions):
        """
        Compute the score according to the heuristic.

        Args:
            predictions (ndarray): Array of predictions

        Returns:
            Array of scores.
        """
        raise NotImplementedError

    def get_uncertainties_generator(self, predictions):
        """
        Compute the score according to the heuristic.

        Args:
            predictions (Iterable): Generator of predictions

        Raises:
            ValueError if the generator is empty.

        Returns:
            Array of scores.
        """
        acc = []
        for pred in predictions:
            acc.append(self.get_uncertainties(pred))
        if len(acc) == 0:
            raise ValueError("No prediction! Cannot order the values!")
        return torch.cat(acc)

    def get_uncertainties(self, predictions):
        """
        Get the uncertainties.

        Args:
            predictions (ndarray): Array of predictions

        Returns:
            Array of uncertainties

        """
        if not isinstance(predictions, Tensor):
            raise TypeError('should be Tensor')

        scores = self.compute_score(predictions)
        scores = self.reduction(scores)
        if not torch.all(torch.isfinite(scores)):
            fixed = 0.0 if self.reversed else 10000
            warnings.warn(f"Invalid value in the score, will be put to {fixed}", UserWarning)
            scores[~torch.isfinite(scores)] = fixed
        return scores

    def reorder_indices(self, scores):
        """
        Order indices given their uncertainty score.

        Args:
            scores (ndarray/ List[ndarray]): Array of uncertainties or
                list of arrays.

        Returns:
            ordered index according to the uncertainty (highest to lowes).

        Raises:
            ValueError if `scores` is not uni-dimensional.
        """
        if isinstance(scores, Sequence):
            scores = torch.cat(scores)

        if scores.ndim > 1:
            raise ValueError(
                (
                    f"Can't order sequence with more than 1 dimension."
                    f"Currently {scores.ndim} dimensions."
                    f"Is the heuristic reduction method set: {self._reduction_name}"
                )
            )
        assert scores.ndim == 1  # We want the uncertainty value per sample.
        ranks = torch.argsort(scores)
        if self.reversed:
            ranks = ranks[::-1]
        ranks = _shuffle_subset(ranks, self.shuffle_prop)
        return ranks

    def get_ranks(self, predictions):
        """
        Rank the predictions according to their uncertainties.

        Args:
            predictions (ndarray): [batch_size, C, ..., Iterations]

        Returns:
            Ranked index according to the uncertainty (highest to lowes).
            Scores for all predictions.

        """
        if isinstance(predictions, types.GeneratorType):
            scores = self.get_uncertainties_generator(predictions)
        else:
            scores = self.get_uncertainties(predictions)

        return self.reorder_indices(scores), scores

    def __call__(self, predictions):
        """Rank the predictions according to their uncertainties.

        Only return the scores and not the associated uncertainties.
        """
        return self.get_ranks(predictions)[0]


class ForceVariance(MyAbstractHeuristic):
    """
    Sort by the highest variance.

    Args:
        shuffle_prop (float): Amount of noise to put in the ranking. Helps with selection bias
            (default: 0.0).
        reduction (Union[str, callable]): function that aggregates the results (default: `mean`).
    """
    def __init__(self, shuffle_prop=0.0, reduction="max"):
        _help = "Need to reduce the output from [n_sample, n_class] to [n_sample]"
        assert reduction != "none", _help
        super().__init__(shuffle_prop=shuffle_prop, reverse=True, reduction=reduction)

    def compute_score(self, predictions):
        assert predictions.ndim >= 3
        predictions = predictions.reshape([predictions.shape[0],-1,3,predictions.shape[-1]])
        pred_mean = predictions.mean(-1)
        devi = torch.sqrt((torch.norm(predictions - pred_mean[...,None], dim=-2)**2).mean(-1)+1e-8)

        return devi
        # return devi.max(-1).values