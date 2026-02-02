#!/usr/bin/env python3
"""
Algorithm 1 (Prange): Attack Simulation for McEliece and HQC Parameter Sets

Author: Anonymous
Date: 2025-11-13
License: MIT
"""

import sys
import math
import numpy as np
from scipy.stats import binom
from helpers import *


class ShiftedBinomial:
    """
    Shifted Binomial distribution model used to simulate weighted probabilities.

    Parameters
    ----------
    n : int
        Number of trials in the binomial distribution.
    p : float
        Probability of success for each trial.
    shift : int
        Constant value added to all samples.

    Attributes
    ----------
    distribution : scipy.stats._distn_infrastructure.rv_frozen
        The underlying binomial distribution object.
    """

    def __init__(self, n: int, p: float, shift: int):
        self.n = n
        self.p = p
        self.shift = shift
        self.distribution = binom(n, p)

    def pmf(self, value: int) -> float:
        """Probability mass function with shift applied."""
        shifted_value = value - self.shift
        return self.distribution.pmf(shifted_value) if shifted_value >= 0 else 0.0

    def cdf(self, value: int) -> float:
        """Cumulative distribution function with shift applied."""
        shifted_value = value - self.shift
        return self.distribution.cdf(shifted_value) if shifted_value >= 0 else 0.0

    def rvs(self, size: int = 1) -> np.ndarray:
        """Random variates from the shifted binomial distribution."""
        samples = self.distribution.rvs(size=size)
        return samples + self.shift


def calc_prob(score: float, d0: ShiftedBinomial, d1: ShiftedBinomial, w: int, n: int) -> float:
    """
    Calculate posterior probability of a given score using weighted binomial models.

    Parameters
    ----------
    score : float
        Observed score value.
    d0 : ShiftedBinomial
        Distribution assuming hypothesis H0.
    d1 : ShiftedBinomial
        Distribution assuming hypothesis H1.
    w : int
        Weight of the error vector.
    n : int
        Code length.

    Returns
    -------
    float
        Posterior probability for the given score.
    """
    numerator = d1.pmf(score) * w
    denominator = d1.pmf(score) * w + d0.pmf(score) * (n - w)
    return numerator / denominator if denominator != 0 else 0.0


def attack(n: int, k: int, w: int, e: np.ndarray, normed: np.ndarray) -> int:
    """
    Simulate a probabilistic attack based on normalized probabilities.

    Parameters
    ----------
    n : int
        Code length.
    k : int
        Dimension of the code.
    w : int
        Weight of the error vector.
    e : np.ndarray
        Error vector.
    normed : np.ndarray
        Normalized probability distribution over indices.

    Returns
    -------
    int
        Sum of recovered error positions.
    """
    indices = np.random.choice(n, n - k, replace=False, p=normed)
    recovered = [e[i] for i in indices]
    return np.sum(recovered)



def main():
    """Main entry point for the simulation."""
    if len(sys.argv) < 4:
        sys.exit("Usage: ./alg1_prange.py <ParameterSet> <NumHints> <NumRuns>")

    parameter = sys.argv[1]
    m = int(sys.argv[2])
    runs = int(sys.argv[3])

    print(f"Simulation of Algorithm 1 (Prange) for parameter set: {parameter}, number of hints: {m}, runs: {runs}")

    n, k, w = get_parameters(parameter)
    total_cost = 0

    for run in range(runs):
        print(f"Run number {run+1}/{runs}")
        # setup instance

        e = error_vec(n, w, parameter)
        H = random_parity_check(n, k)
        # compute perfect hints
        z = H @ e
        # take only m hints
        H_m = H[:m]
        z_m = z[:m]

        psis = [(i, psi(i, H_m, z_m, w, m)) for i in range(n)]

        # Define shifted binomial models
        p = 0.5
        d0 = ShiftedBinomial(m * w, p, 0)
        d1 = ShiftedBinomial(m * (w - 1), p, m)

        # Compute normalized probabilities
        probs = [(x[0], calc_prob(x[1], d0, d1, w, n)) for x in psis]
        total = np.sum([x[1] for x in probs])
        normed = np.array([x[1] / total for x in probs])

        # Repeated attack trials
        best, count = 0, 0
        while best < w:
            count += 1
            res = attack(n, k, w, e, normed)
            if res > best:
                best = res
                print(f"Best try has {w - best} missed error coordinates - 2^{log2(count)} iterations.")

        print(f"({m}, 2^{math.log2(count):.2f})")
        total_cost += count

    # Save average results
    avg_cost = total_cost / runs
    with open(f"alg1-prange-{parameter}.txt", "a") as file:
        file.write(f"({m}, 2^{math.log2(avg_cost):.2f})\n")


if __name__ == "__main__":
    main()
