#!/usr/bin/env python3
"""
Algorithm 2 (Prange): Derivation of optimal alpha (number of fixed positions) by computing the costs
and taking the alpha with minimal cost for each parameter setting and hint number.

Author: Anonymous
Date: 2025-11-13
License: MIT
"""

import sys
import pprint
from helpers import *
from math import comb, log2
import pprint

if __name__ == "__main__":
    
    # read command-line args or set defaults
    parameter = str(sys.argv[1]) if len(sys.argv) > 1 else "McEliece1"
    form = str(sys.argv[2]) if len(sys.argv) > 2 else "systematic"
    alg = "prange"

    print("Computing the optimal number of fixed positions, the expected number of ones and the corresponding cost for Alg2.")
    print(f"parameter: {parameter}, matrix: {form}, algorithm: {alg}.")

    # get code parameters
    n, k, w = get_parameters(parameter)

    # setup instance
    e = error_vec(n, w, parameter)
    H = parity_check(n, k, form)
    
    # evaluation steps for m
    if parameter in HQC:
        steps = list(range(0, 7000, 100))
    else:
        steps = list(range(0, n-k+1, 50)) + [n-k]


    best_results = {}

    # alpha resolution differs by scheme (for efficiency)
    if parameter in HQC:
        alpha_steps = list(range(0, n-k+1, 50)) + [n-k]
    elif parameter in MCELIECE:
        alpha_steps = list(range(0, n-k+1, 10)) + [n-k]

    # loop over number of fixed positions
    for m in steps:
        best_cost = float("inf")
        best_alpha = None

        # try all alpha values
        for alpha in alpha_steps:
            
            # random parity-check matrix path
            if form == "random":
                eo = int(exp_num_ones(m, alpha, n, w))
                cost = prange_cost(alpha, n, k, w, eo)
                #print(f"alpha:{alpha}, eo:{eo}, cost:2^{cost} ")

            # systematic parity-check matrix
            elif form == "systematic":
                dim = k                       # size of random block
                weight = round(w * k / n)     # expected weight in random block
                eo = int(exp_num_ones(m, alpha, dim, weight)) # expected number of ones in random block
                cost, perm_l, perm_r, w_l, w_r_perm = systematic_prange_cost(alpha, n, k, w, eo, parameter)
                #print(f"cost:{cost}, alpha:{alpha}, eo:{eo}, left:{perm_l}, right:{perm_r}, w_l:{w_l}, w_r:{w_r_perm}")
            # track minimum cost
            if cost < best_cost:
                best_cost = cost
                best_alpha = alpha
                best_eo = eo
                if form == "systematic":
                    best_perm_l = perm_l
                    best_perm_r = perm_r
                    best_w_l = w_l
                    best_w_r_perm = w_r_perm

        # store optimal values for this m
        if form == "systematic":
            best_results[m] = (best_alpha, best_eo, best_cost, best_perm_l, best_perm_r, best_w_l, best_w_r_perm)
        else:
            best_results[m] = (best_alpha, best_eo, best_cost)

        print(f"({m},{best_results[m][2]})")

    # append results to EXPECTED_PRANGE.py
    with open("EXPECTED_PRANGE.py", "a") as f:
        f.write(f"{parameter}_{form}_{alg} =")
        pprint.pprint(best_results, stream=f)
        f.write("\n")
