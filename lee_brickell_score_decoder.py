#!/usr/bin/env python3
"""
Algorithm 2 (Stern): Attack Simulation for McEliece and HQC Parameter Sets

Author: Anonymous
Date: 2025-11-13
License: MIT
"""

import numpy as np
import sys
from helpers import *
from EXPECTED_STERN import *

P_LEE_BRICKELL = 3

def lee_brickell_score(H, z, e, w, m, n ,k):
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

    # if solvable instance return True
    if ones_initial >= w - P_LEE_BRICKELL:
        return True
    # else return False (not solvable)            
    return False

def bin_noise(d):
    noise = -d + np.random.binomial(2*d, 0.5)
    return noise

if __name__ == "__main__":
    parameter = str(sys.argv[1]) if len(sys.argv) > 1 else "McEliece1"
    form = str(sys.argv[2]) if len(sys.argv) > 2 else "systematic"
    noise = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    runs = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    alg = "stern"

    print("Re-Implementation of Lee-Brickell Score Decoder Attack. Outputs the fraction of successful runs for a given number of hints m.")
    print(f"parameter: {parameter}, matrix: {form}, algorithm: {alg}, runs:{runs}")

    n, k, w = get_parameters(parameter)

    # set stepsize dependent on the scheme
    if parameter in HQC:
        steps = list(range(0,n-k+1,100)) + [n-k]
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
            # add noise
            noise_vec = [bin_noise(noise) for _ in range(n-k)]
            z += noise_vec
            # take only m hints
            H_m = H[:m]
            z_m = z[:m]
            # run the algorithm
            success = lee_brickell_score(H_m, z_m, e, w, m, n, k)
            # write to stdout
            print(f"m: {m}, Attack applicable: {success}")
            # write to results
            if m not in results:
                results[m] = 1 if success else 0
            else:
                results[m] += 1 if success else 0

        success_rate = results[m] / runs

        with open(f"lee-brickell-scorer-successrate-{parameter}-{form}-{noise}.txt", "a") as file:
            file.write(f"({m}, {success_rate})\n")





