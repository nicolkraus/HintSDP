# Syndrome Decoding with Perfect Hints

This repository contains the source code for Algorithm 1, Algorithm 2, a script to compute the expected cost of Algorithm 2 as well as a script to run the (Integer) Linear Program and the Lee-Brickell score decoder.

## Parameters

The scripts take the following command-line parameters.

#### Parameter set

One of the predefined code parameter sets:

- `McEliece1`
- `McEliece2`
- `McEliece3`
- `McEliece4`
- `McEliece5`
- `HQC1`
- `HQC2`
- `HQC3`

Each identifier corresponds to a specific \((n, k, w)\) triple, where `n` is the code length, `k` the code dimension and `w` the error weight as defined in the paper.

> Note: Computations for HQC parameter sets require significant resources and are best run on a server.

#### Matrix form

Determines how the parity-check matrix is generated:

- `random` – random parity-check matrix  
- `systematic` – systematic parity-check matrix

#### Number of hints m

Specifies how many perfect hints are provided.  Must be less than or equal to `n-k`.

#### Noise parameter d

Specifies the noise level. For each hint, a noise term e ~ Bin(2*d, 0.5)-d is sampled and added.

#### Runs

Defines how many independent executions are performed. Each run uses a fresh error vector and parity-check matrix.


## Algorithm 1

To reproduce the results of Algorithm 1, run:

```bash
python3 alg1_stern.py [parameter] [number_of_hints_m] [noise d] [runs]
```


## Algorithm 2

### Expected cost computation

To compute expected costs and optimal choices of the number of fixed positions:

```bash
python3 get_expected_stern.py [parameter] [matrix_form] [noise d]
```

> Note: We precomputed the expected costs and optimal number of fixed positions, so the above script does not need to be run for the experiments.

### Experimental evaluation

To run the experiments for Algorithm 2:

```bash
python3 alg2_stern.py [parameter] [matrix_form] [noise d] [runs]
```

## Prange with Hints
To run the experiments for Prange's algorithm with hints, change into directory `prange` and run the experiments in the same way as above for Stern:

```bash
python3 alg1_prange.py [parameter] [number_of_hints_m] [runs]
python3 alg2_prange.py [parameter] [matrix_form] [runs]
```


## Integer Linear Programming (ILP)

To reproduce the ILP results from the submission, run:

```bash
python3 ILP.py [parameter] [matrix_form]
```
> Note: This script requires an active Gurobi license. For academic use, the license can be obtained for free.
