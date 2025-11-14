import numpy as np
import random
from math import log2, comb
from scipy.stats import binom

# Constant representing the cost threshold for the Stern algorithm.
# For all supported parameter sets, the cost of comb(k/2, P_STERN/2)
# is smaller than the cost of Gaussian elimination.
P_STERN = 6

# Supported parameter sets for HQC and McEliece schemes.
HQC = ["HQC1", "HQC2", "HQC3"]
MCELIECE = ["McEliece1", "McEliece2", "McEliece3", "McEliece4", "McEliece5"]


def weighted_vec(n, w):
    """
    Generate a binary vector of length n with exactly w entries set to 1.

    Parameters
    ----------
    n : int
        Length of the vector.
    w : int
        Hamming weight (number of 1s) in the vector.

    Returns
    -------
    numpy.ndarray
        Binary vector of length n and weight w.
    """
    vector = np.zeros(n, dtype=int)  # Initialize all entries to zero
    indices = np.random.choice(n, size=w, replace=False)  # Randomly select w unique positions
    vector[indices] = 1  # Set those positions to 1
    return vector


def weighted_HQC_vec(n, w):
    """
    Generate a binary HQC-style vector of length n with weight w,
    where half the weight is placed in each half of the vector.

    Parameters
    ----------
    n : int
        Length of the vector (even number).
    w : int
        Total Hamming weight.

    Returns
    -------
    numpy.ndarray
        Binary vector of length n and total weight w, split evenly between halves.
    """
    vector = np.zeros(n, dtype=int)  # Initialize vector of zeros
    # Select w/2 indices from the first half
    indices1 = np.random.choice(n // 2, size=w // 2, replace=False)
    # Select w/2 indices from the second half (offset by n/2)
    indices2 = [x + n // 2 for x in np.random.choice(n // 2, size=w // 2, replace=False)]
    # Set both selected halves to 1
    vector[indices1] = 1
    vector[indices2] = 1
    return vector


def binary_uniform_vec(n):
    """
    Generate a uniformly random binary vector of length n.

    Parameters
    ----------
    n : int
        Length of the binary vector.

    Returns
    -------
    numpy.ndarray
        Binary vector with entries uniformly chosen from {0, 1}.
    """
    v = np.array([random.randrange(0, 2) for _ in range(n)])  # Independent uniform bits
    return v


def systematic_parity_check(n, k):
    """
    Generate a systematic-form parity-check matrix of size (n-k) x n.

    H = [I | R], where I is the identity matrix of size (n-k)
    and R is a random binary matrix of size (n-k) x k.

    Parameters
    ----------
    n : int
        Code length.
    k : int
        Code dimension.

    Returns
    -------
    numpy.ndarray
        Systematic parity-check matrix.
    """
    I = np.eye(n - k, dtype=int)  # Identity matrix for the systematic part
    R = np.random.randint(0, 2, size=(n - k, k))  # Random binary matrix
    H = np.concatenate((I, R), axis=1)  # Concatenate horizontally
    return H


def random_parity_check(n, k):
    """
    Generate a fully random parity-check matrix of size (n-k) x n.

    Parameters
    ----------
    n : int
        Code length.
    k : int
        Code dimension.

    Returns
    -------
    numpy.ndarray
        Random binary parity-check matrix.
    """
    H = np.random.randint(0, 2, size=(n - k, n))  # Random 0/1 entries
    return H


def psi(index, H, z, w, m):
    """
    Compute the score ψ(i) for a given index in the attack algorithm.
    The score measures the correlation between a column of H and the vector z.

    Parameters
    ----------
    index : int
        Index of the column in H to evaluate.
    H : list[numpy.ndarray]
        List of binary sample vectors.
    z : list[float]
        Corresponding inner product values (sample results).
    w : int
        Weight parameter of the target error vector.
    m : int
        Number of samples.

    Returns
    -------
    float
        The computed ψ-score for column `index`.
    """
    col = [h[index] for h in H]  # Extract the column values for the given index
    compl_col = np.ones(m) - col  # Complement of the column (flip bits)
    compl_int_syndr = w * np.ones(m) - z  # Complement of the syndrome vector
    score = np.dot(z, col) + np.dot(compl_int_syndr, compl_col)
    return score


def prange_cost(alpha, n, k, w, eo):
    """
    Compute the binary logarithm of the expected cost of a Prange decoding step.

    Parameters
    ----------
    alpha : int
        Number of fixed positions.
    n : int
        Code length.
    k : int
        Code dimension.
    w : int
        Weight of the target error vector.
    eo : int
        Expected number of errors in the alpha fixed positions.

    Returns
    -------
    float
        Log2 of the expected number of operations, or infinity if infeasible.
    """
    # Boundary cases: trivial or impossible decoding
    if w - eo < 1:
        cost = 0
    elif n - k - alpha < w - eo:
        cost = float("inf")
    else:
        # Cost formula derived from Prange’s combinatorial decoding model
        cost = log2(comb(n - alpha, w - eo)) - log2(comb(n - k - alpha, w - eo))
    return cost

def exp_num_ones(m, alpha, n, w):
    """
    Estimate the expected number of ones in highest scored alpha positions.

    Parameters
    ----------
    m : int
        Number of samples (measurements or hints) available.
    alpha : int
        Target number of bits (or positions) satisfying a certain condition.
    n : int
        Code length.
    w : int
        Hamming weight of the error vector.

    Returns
    -------
    float
        Estimated expected number of ones (E[#1s]) for the given parameters.
    """
    if m <= 0:
        return 0

    # Define function f(B) whose root determines the correct threshold B
    f = lambda B: alpha - (n - w) * (1 - binom.cdf(B, m * w, 1 / 2)) \
                  - w * (1 - binom.cdf(B - m, m * (w - 1), 1 / 2))

    # Binary search bounds for B
    lo, hi = 0, m * w + m

    # Perform binary search for 60 iterations to approximate the root
    for _ in range(60):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid

    # Midpoint is our best approximation of B
    B = (lo + hi) / 2

    # Expected number of ones from the binomial tail probability
    exp_num_ones = w * (1 - binom.cdf(B - m, m * (w - 1), 1 / 2))
    return exp_num_ones


def stern_cost(alpha, n, k, w, eo):
    """
    Compute the expected binary logarithmic cost of Stern's decoding algorithm.

    Parameters
    ----------
    alpha : int
        Reduction parameter (defines subset of positions considered).
    n : int
        Code length.
    k : int
        Code dimension.
    w : int
        Weight of the error vector.
    eo : int
        Number of errors already obtained (partial recovery).

    Returns
    -------
    float
        Log2 of the expected cost for the Stern algorithm, or infinity if
        parameters make the decoding infeasible.
    """
    if w - eo - P_STERN < 1:
        cost = 0 # corresponding to 2^0
    elif n - k - alpha < w - eo - P_STERN:
        cost = float("inf")
    else:
        # Stern's cost formula (logarithmic form)
        cost = log2(comb(n - alpha, w - eo)) \
             - log2(comb(n - k - alpha, w - eo - P_STERN)) \
             - log2(comb(k, P_STERN))
    return cost



def get_parameters(parameter_name: str):
    """
    Retrieve (n, k, w) parameters for supported McEliece and HQC schemes.

    Parameters
    ----------
    parameter_name : str
        Name of the parameter set (e.g., "McEliece1", "HQC2").

    Returns
    -------
    tuple[int, int, int]
        (n, k, w) parameter triple.

    Raises
    ------
    ValueError
        If the given parameter set name is not recognized.
    """
    params = {
        "McEliece1": (3488, 2720, 64),
        "McEliece2": (4608, 3360, 96),
        "McEliece3": (6688, 5024, 128),
        "McEliece4": (6960, 5413, 119),
        "McEliece5": (8192, 6528, 128),
        "HQC1": (35398, 17699, 132),
        "HQC2": (71702, 35851, 200),
        "HQC3": (115274, 57637, 262),
    }

    if parameter_name not in params:
        raise ValueError("Invalid parameter set. Choose from McEliece[1–5] or HQC[1–3].")

    return params[parameter_name]


def error_vec(n, w, parameter):
    """
    Generate an error vector suitable for the given parameter set.

    Parameters
    ----------
    n : int
        Length of the error vector.
    w : int
        Hamming weight.
    parameter : str
        Parameter set name (determines HQC vs McEliece style).

    Returns
    -------
    numpy.ndarray
        Error vector of weight w appropriate to the scheme.
    """
    if parameter in HQC:
        return weighted_HQC_vec(n, w)
    else:
        return weighted_vec(n, w)


def parity_check(n, k, form):
    """
    Generate a parity-check matrix in either random or systematic form.

    Parameters
    ----------
    n : int
        Code length.
    k : int
        Code dimension.
    form : str
        Either "random" or "systematic" — determines matrix structure.

    Returns
    -------
    numpy.ndarray
        Parity-check matrix.
    """
    if form == "random":
        return random_parity_check(n, k)
    elif form == "systematic":
        return systematic_parity_check(n, k)
    
# Computes the cost of Stern’s algorithm in the systematic setting.
# Columns are chosen proportionally from the random part (size k) and the identity part (size n–k),
# based on the expected distribution of error positions.
def systematic_stern_cost(alpha, n, k, w, eo):
    w_new = w - eo # Total remaining error weight after enumerating eo coordinates
    w_r = round(w * k / n) # Expected number of error positions in the random part (size k)
    w_l = w - w_r # Expected number of error positions in the identity part (size n-k)
    w_r_perm = w_r - eo     # Expected remaining error positions (after eo already handled) in random part
    share_r = w_r_perm / w_new # Fraction of remaining weight that lies in the random part
    perm_r = round((n - k - alpha) * share_r) # Number of columns to draw from the random part (k-alpha remaining)
    perm_l = (n - k - alpha) - perm_r # Remaining columns must come from the identity part
    get_r = w_r_perm - P_STERN # Required number of errors we must obtain from the chosen coordinates of the random part

    # Case 1: Less than P_STERN errors in the random part
    if get_r < 0 and (n - k - alpha) >= w_l:
        cost = (
            log2(comb(n - k, w_l))      # ways to choose identity-part errors
            - log2(comb(n - k - alpha, w_l)) # remove alpha excluded coordinates
        )

    # Case 2: Infeasible configurations
    elif (
        (n - k - alpha) < (w - eo - P_STERN)    # too few remaining coordinates
        or (w_r_perm > perm_r and get_r >= 0)   # insufficient random-part columns to extract required errors
        or (w_l > perm_l)                  # insufficient identity-part columns to extract required errors
    ):
        cost = float("inf")

    # Case 3: Standard feasible cost computation
    else:
        cost = (
            log2(comb(n - k, w_l))               # all possible error positions in identity part
            + log2(comb(k - alpha, w_r_perm))    # all possible error positions in random part
            - log2(comb(perm_r, get_r))          # good choices (random part)
            - log2(comb(perm_l, w_l))            # good choices (identity part)
            - log2(comb(k, P_STERN))             # enumeration of P_STERN error positions
        )

    return cost, perm_l, perm_r

    
# Computes the cost of Stern’s algorithm in the systematic setting when all parameters are given.
# Allows to compute the expected cost for single experiments
def experiments_systematic_stern_cost(alpha, n, k, w, fixed_ones, perm_l, perm_r, w_l, w_r_perm):
    get_r = w_r_perm - P_STERN # Required number of errors we must obtain from the chosen coordinates of the random part
    
    # Case 1: Less than P_STERN errors in the random part
    if get_r < 0 and (n - k - alpha) >= w_l:
        cost = (
            log2(comb(n - k, w_l))      # ways to choose identity-part errors
            - log2(comb(n - k - alpha, w_l)) # remove alpha excluded coordinates
        )

    # Case 2: Infeasible configurations
    elif (
        (n - k - alpha) < (w - fixed_ones - P_STERN)    # too few remaining coordinates
        or (w_r_perm > perm_r and get_r >= 0)   # insufficient random-part columns to extract required errors
        or (w_l > perm_l)                  # insufficient identity-part columns to extract required errors
    ):
        print(w_r_perm, perm_r )
        cost = float("inf")

    # Case 3: Standard feasible cost computation
    else:
        cost = (
            log2(comb(n - k, w_l))               # all possible error positions in identity part
            + log2(comb(k - alpha, w_r_perm))    # all possible error positions in random part
            - log2(comb(perm_r, get_r))          # good choices (random part)
            - log2(comb(perm_l, w_l))            # good choices (identity part)
            - log2(comb(k, P_STERN))             # enumeration of P_STERN error positions
        )

    return cost

