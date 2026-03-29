import numpy as np
import warnings
from typing import Tuple

def duplicated_tol(x, tol=1e-10, from_last=False):
    """
    Tolerance-aware version of R's duplicated().

    Parameters
    ----------
    x : array-like
        Input array.
    tol : float
        Absolute tolerance for equality.
    from_last : bool
        If True, mark duplicates keeping the *last* occurrence
        (equivalent to duplicated(..., fromLast = TRUE) in R).

    Returns
    -------
    dup : ndarray of bool
        True where values are duplicates.
    """
    x = np.asarray(x)
    n = len(x)
    dup = np.zeros(n, dtype=bool)

    seen = []

    if from_last:
        it = range(n - 1, -1, -1)
    else:
        it = range(n)

    for i in it:
        if any(abs(x[i] - s) < tol for s in seen):
            dup[i] = True
        else:
            seen.append(x[i])

    return dup

def compute_lqd_cut(
        lqd_curve: np.array,
        k: int = 15,
        threshold: float = 7.0,
        max_frac: float = 0.1,
        cut_nan: bool = False
    ) -> Tuple[int, int]:
    """
    Compute adaptive boundary cuts for a single LQD curve that has problematic values (extreme or nan).

    Parameters
    ----------
    lqd_curve : array-like, shape (M,)
        One LQD realization.
    k : int
        Number of boundary points to inspect on each side.
    threshold : float
        LQD value above which points are considered unstable.
    cut_nan : bool
        Include cutting NaN values that appear in the left or right endpoints of the support

    Returns
    -------
    (left_cut, right_cut) : tuple of ints
    """
    arr = np.asarray(lqd_curve)
    M = arr.shape[0]

    nan_left = 0
    nan_right = 0

    # ------------------------------
    # Optional NaN trimming
    # ------------------------------
    if cut_nan:
        finite = np.isfinite(arr)

        if not finite.any():
            raise ValueError("All LQD values are invalid.")

        nan_left = int(np.argmax(finite))
        nan_right = int(np.argmax(finite[::-1]))

        if nan_left > 0 or nan_right > 0:
            warnings.warn(
                f"LQD contained NaN/inf at boundaries — cutting "
                f"{nan_left} left and {nan_right} right points.",
                RuntimeWarning,
                stacklevel=2
            )

    # ------------------------------
    # Threshold trimming AFTER NaNs
    # ------------------------------
    left_slice = arr[nan_left:nan_left + k]
    right_slice = arr[M - nan_right - k:M - nan_right]

    left_thr = int(np.sum(left_slice > threshold))
    right_thr = int(np.sum(right_slice > threshold))

    left = nan_left + left_thr
    right = nan_right + right_thr

    # Safety cap
    max_cut = int(max_frac * M)

    return min(left, max_cut), min(right, max_cut)

def duplicated_tol_sorted(x, tol=1e-10, from_last=False):
    """
    Tolerance-aware version of R's duplicated() for sorted vectors.

    Parameters
    ----------
    x : array-like
        Input array (assumed monotone).
    tol : float
        Absolute tolerance for equality.
    from_last : bool
        If True, mark duplicates keeping the *last* occurrence.

    Returns
    -------
    dup : ndarray of bool
        True where values are duplicates.
    """
    x = np.asarray(x)

    if len(x) == 0:
        return np.zeros(0, dtype=bool)

    if from_last:
        dx = np.abs(np.diff(x[::-1])) < tol
        return np.r_[dx, False][::-1]
    else:
        dx = np.abs(np.diff(x)) < tol
        return np.r_[False, dx]
