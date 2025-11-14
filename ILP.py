#!/usr/bin/env python3
"""
Integer Linear Programming: Attack Simulation for McEliece and HQC Parameter Sets

Author: Nicolai Kraus
Date: 2025-11-13
License: MIT
"""

import numpy as np
import sys
from gurobipy import GRB, quicksum, Model
from helpers import *

# Solve the ILP for given H, z, e, w
def ilp(H, z, e, w, parameter, n):
    correct = 0
    # Setup model
    model = Model()
    model.setParam('OutputFlag', 0)
    # Variables
    x_vars = []
    for i in range(n):
        x = model.addVar(lb=0 , ub=1 , vtype=GRB.CONTINUOUS, name=f'x_{i}')
        x_vars.append(x)
    # Constraints
    for sig, coeffs in zip(z, H):
        expr = quicksum(coeff * var for coeff, var in zip(coeffs, x_vars))
        model.addConstr(expr == sig)
    # Additional weight constraint for HQC
    if parameter in HQC:
        model.addConstr(quicksum(x_vars[i] for i in range(n//2)) == w//2, "first_half_ones")
        model.addConstr(quicksum(x_vars[i] for i in range(n//2, n)) == w//2, "second_half_ones")
    # weight constraint for McEliece
    if parameter in MCELIECE:
        model.addConstr(quicksum(x_vars[i] for i in range(n)) == w, "total_ones")
    # Objective (minimize 0)
    model.setObjective(quicksum(0 * var for var in x_vars), GRB.MINIMIZE)
    # Solve
    model.optimize()
    # Extract the estimated b vector from the model coefficients
    if model.Status == GRB.OPTIMAL:
        e_est = np.array([var.X for var in x_vars])
        correct = np.sum(np.round(e_est) == e)

    return correct

def main():
    """Main entry point for the simulation."""
    parameter = str(sys.argv[1]) if len(sys.argv) > 1 else "HQC1"
    form = str(sys.argv[2]) if len(sys.argv) > 2 else "systematic"
    stepsize = int(sys.argv[3]) if len(sys.argv) > 3 else 100 

    print(f"parameter: {parameter}, matrix is {form}, stepsize: {stepsize}")

    n, k, w = get_parameters(parameter)

    e = error_vec(n, w, parameter)        
    H = parity_check(n, k, form)

    # compute perfect hints
    z = H @ e

    # list with steps of stepsize and the full integer syndrome
    steps = list(range(stepsize,n-k+1,stepsize))+[n-k]
    # increase number of hints until secret is found
    for m in steps:
        # take only m hints
        H_m = H[:m]
        z_m = z[:m]
        # run the algorithm
        res = ilp(H_m, z_m, e, w, parameter, n)
        # write to stdout
        print(f"m={m}, correct positions={res}/{n}")
        # write results to file
        with open(f"ILP-{parameter}-{form}-{stepsize}.txt", "a") as file:
            file.write(f"({m}, {res}/{n})\n")
        # stop if secret is found
        if res == n:
            break

if __name__ == "__main__":
    main()
   




