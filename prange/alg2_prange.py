#!/usr/bin/env python3
"""
Algorithm 2 (Prange): Attack Simulation for McEliece and HQC Parameter Sets

Author: Anonymous
Date: 2025-11-13
License: MIT
"""

import numpy as np
import sys
from helpers import *
from EXPECTED_PRANGE import *

def alg2(H, z, e, w, m, n ,k, form):
    # select algorithm dependent on the form of the parity-check matrix
    if form == "systematic":
        return alg2_systematic(H, z, e, w, m, n, k)
    else:
        return alg2_random(H, z, e, w, m, n, k)

def alg2_random(H, z, e, w, m, n ,k):
    # access precomputed dictionary dependent on parameter
    dict = f"{parameter}_random_{alg}"
    # fetch optimal alpha, expected ones and cost
    alpha, _, _ = globals()[dict][m]

    # score the indices
    scores = []
    for i in range(n):
        t = (i, psi(i,H,z,w,m))
        scores.append(t)
    # sort scores in descending order
    scores.sort(key = lambda x:x[1], reverse = True)   

    # check if problem is trivially solved in poly time
    ids_initial = [x[0] for x in scores][:n-k]
    ones_initial = int(np.sum([e[i] for i in ids_initial]))
    if ones_initial == w:
        return 0 # corresponding to 2^0
    
    # else fix alpha top positions
    ids = [x[0] for x in scores][:alpha]
    ones = int(np.sum([e[i] for i in ids]))
    cost = prange_cost(alpha,n,k,w,ones)
                
    return cost


def alg2_systematic(H, z, e, w, m, n ,k):
    # access precomputed dictionary dependent on parameter
    dict = f"{parameter}_systematic_{alg}"
    # fetch optimal alpha, expected ones and cost
    alpha, _, _, perm_l, perm_r, _, _ = globals()[dict][m]
    # expected weight in the random part
    rand_weight = round(w*k/n)

    #score random part
    scores = []
    for i in range(n-k,n):
        t = (i, psi(i,H,z,rand_weight,m))
        scores.append(t)
    # sort scores
    scores.sort(key = lambda x:x[1], reverse = True)   

    # check if problem is trivially solved in poly time
    ids_initial = [x[0] for x in scores][:w//2]
    ones_initial = int(np.sum([e[i] for i in ids_initial]))
    if ones_initial == w//2: # all coordinates corresponding to the random part are found
        return 0 # corresponding to 2^0
    
    # fix alpha top positions
    id_random = [x[0] for x in scores][:alpha]
    # get fixed and free weights dependent on the number of alpha fixed positions (to compute the remaining expected cost)
    fixed_ones = int(np.sum([e[i] for i in id_random]))
    w_r = int(np.sum(e[i] for i in range(n-k,n))) - fixed_ones
    w_l = int(np.sum(e[i] for i in range(n-k)))
    #print(f"w_r:{w_r}, w_l:{w_l}, fixed_ones:{fixed_ones}, alpha:{alpha}")

    cost = experiments_systematic_prange_cost(alpha, n, k, w, fixed_ones, perm_l, perm_r, w_l, w_r, parameter)

    return cost

if __name__ == "__main__":
    parameter = str(sys.argv[1]) if len(sys.argv) > 1 else "McEliece1"
    form = str(sys.argv[2]) if len(sys.argv) > 2 else "random"
    runs = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    alg = "prange"

    print(f"parameter: {parameter}, matrix: {form}, algorithm: {alg}, runs:{runs}")

    n, k, w = get_parameters(parameter)

    # set stepsize dependent on the scheme
    if parameter in HQC:
        steps = list(range(0,7000,100))
    elif parameter in MCELIECE:
        steps = list(range(0,n-k+1,50)) + [n-k]

    # setup results dictionary
    results = {}
    best80percent = {}

    # increase number of hints until secret is found
    for m in steps:
        # run multiple times and average results
        for run in range(runs):
            print(f"Run {run+1}/{runs}")
            # setup instance
            e = error_vec(n, w, parameter)
            H = parity_check(n, k, form)
            # compute perfect hints
            z = H @ e
            # take only m hints
            H_m = H[:m]
            z_m = z[:m]
            # run the algorithm
            cost = alg2(H_m, z_m, e, w, m, n, k, form)
            # write to stdout
            print(f"m: {m}, cost: 2^{cost}")
            # write to results
            if m not in results:
                results[m] = []
            results[m].append(cost)

        cost = sum(results[m]) / len(results[m])
        best80percent[m] = sorted(results[m])[round(runs*0.8)-1]

        with open(f"alg2-prange-{parameter}-{form}.txt", "a") as file:
            file.write(f"({m}, 2^{cost})\n")

    # Find the m whose cost is 0 in 80% of the executions and write file
    poly_ms = [m for m, v in best80percent.items() if v == 0]
    with open(f"alg2-prange-{parameter}-{form}.txt", "a") as file:
        file.write(f"Poly time in 80% of the executions for {poly_ms}\n")





