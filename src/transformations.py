import numpy as np
import warnings
from scipy import integrate
from scipy import interpolate 
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.interpolate import interp1d, CubicSpline, InterpolatedUnivariateSpline
from scipy.integrate import cumulative_trapezoid, quad, trapezoid
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from typing import Tuple

# Fixed M to better approximate functions. Smaller M (like 256) generates worse approximations.
# M = 5001

#---------------------------------------------------------------------
#----------------------------- UTILITIES -----------------------------
#---------------------------------------------------------------------
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

#---------------------------------------------------------------------
#----------------------- TRANSFORMATIONS -----------------------------
#---------------------------------------------------------------------

def dens2lqd(dens, dSup, lqdSup=None, t0=None, verbose=True):
    """
    ----------------------------------------------------------------------
    Mathematic Background
    ----------------------------------------------------------------------
    Let f be a density supported on an interval [a, b], with cumulative
    distribution function (CDF):

        F(x) = ∫_a^x f(s) ds.

    Define the quantile function:

        Q(u) = F^{-1}(u),     u ∈ [0, 1].

    The *quantile-density function* is:

        q(u) = Q'(u) = 1 / f(Q(u)).
s
    The *log-quantile-density transform* is

        L(u) = -log f(Q(u)) = log q(u).                    (1)

    As shown in Kokoszka & Reimherr (2019), the LQD map:

        f  ↦  L

    takes densities from a nonlinear manifold into an (approximately) 
    linear Hilbert space L^2[0,1], enabling classical FDA tools.

    Parameters
    ----------
    dens : array_like
        A 1D array of density values. Must be non-negative.  
        If the integral is not 1, the function renormalizes it.

    dSup : array_like
        The support points corresponding to `dens`. Must be strictly increasing.

    lqdSup : array_like, optional
        Desired evaluation grid for the LQD, assumed to lie in `[0,1]`.  
        If not provided, defaults to a uniform grid of the same size as `dSup`.  
        If the grid does not span `[0,1]`, it is replaced by the default.

    t0 : float, optional
        A reference location in the original support where the CDF-based
        constant `c` is evaluated.  
        If not inside `dSup`, it is replaced by the closest support value.  
        Default is the lower boundary `dSup[0]`.

    verbose : bool, default=True
        Whether to print warnings about renormalization, support truncation,
        and boundary corrections.

    Returns
    -------
    lqdSup : ndarray
        The LQD support grid (in `[0,1]`).

    lqd : ndarray
        The log-quantile-density values evaluated on `lqdSup`.

    c : float
        A constant equal to the CDF evaluated at `t0`, i.e.
        `c = ∫_{dSup[0]}^{t0} dens(s) ds`.

    Notes
    -----
    - Zero density values cause the support to be truncated to ensure
      strictly positive densities, as the LQD is undefined at zeros.
    - CDF monotonicity issues (often appearing in KDE-based estimates)
      are corrected by discarding duplicated CDF entries to ensure
      a valid interpolation grid.
    - If the transformed LQD has infinite values at boundaries, interpolation
      is performed excluding those endpoints, which are assigned `+inf`.
    """
    dens = np.asarray(dens)
    dSup = np.asarray(dSup)

    # Default t0
    if t0 is None:
        t0 = dSup[0]

    # ---- Check density requirements ----
    if np.any(dens < 0):
        raise ValueError("Please correct negative density values.")

    if abs(trapezoid(dens, dSup) - 1) > 1e-5:
        if verbose:
            print("Density does not integrate to 1 with tolerance 1e-5 - renormalizing now.")
        dens = dens / trapezoid(dens, dSup)

    # ---- Handle zero density values by truncating support ----
    if np.any(dens == 0):
        if verbose:
            print("There are some zero density values - truncating support grid so all are positive")

        positive_idx = np.where(dens > 0)[0]
        lbd, ubd = positive_idx[0], positive_idx[-1]

        dens = dens[lbd:ubd+1]
        dSup = dSup[lbd:ubd+1]

        dens = dens / trapezoid(dens, dSup)

    N = len(dSup)

    # ---- Check LQD output grid ----
    if lqdSup is None:
        lqdSup = np.linspace(0, 1, N)
    else:
        lqdSup = np.asarray(lqdSup)
        if not (np.isclose(lqdSup.min(), 0) and np.isclose(lqdSup.max(), 1)):
            if verbose:
                print("Problem with support of the LQD domain’s boundaries - resetting to default.")
            lqdSup = np.linspace(0, 1, N)

    # ---- Check t0 ----
    if t0 not in dSup:
        if verbose:
            print("t0 is not a value in dSup - resetting to closest value")
        t0 = dSup[np.argmin(np.abs(dSup - t0))]

    M = len(lqdSup)
    c_ind = np.where(dSup == t0)[0][0]

    # ---- Compute CDF and constant c ----
    tmp = cumulative_trapezoid(dens, dSup, initial=0)
    c = tmp[c_ind]

    # ---- Remove duplicated CDF values (monotonicity issues in KDE) ----
    left_dup  = np.concatenate([np.diff(tmp[:N//2]) == 0, [False]])
    right_dup = np.concatenate([[False], np.diff(tmp[N//2:]) == 0])

    # NOTE: In R: !c(indL, indR)
    keep = ~(np.concatenate([left_dup, right_dup]))

    qtemp = tmp[keep]
    lqd_temp = -np.log(dens[keep])

    # ---- Interpolate lqd on the desired LQD support ----
    lqd = np.zeros(M)

    # Handle infinite boundary values
    temp_first_inf = np.isinf(lqd_temp[0])
    temp_last_inf = np.isinf(lqd_temp[-1])

    if temp_first_inf or temp_last_inf:

        tmpInd = np.arange(len(qtemp))
        Ind = np.arange(M)

        if temp_first_inf:
            lqd[0] = np.inf
            tmpInd = tmpInd[1:]
            Ind = Ind[1:]

        if temp_last_inf:
            lqd[-1] = np.inf
            tmpInd = tmpInd[:-1]
            Ind = Ind[:-1]

        interp = interp1d(qtemp[tmpInd], lqd_temp[tmpInd],
                          kind="linear", fill_value="extrapolate")
        lqd[Ind] = interp(lqdSup[Ind])

    else:
        interp = interp1d(qtemp, lqd_temp, kind="linear",
                          fill_value="extrapolate")
        lqd = interp(lqdSup)

    return lqdSup, lqd, c

def obtain_lqds(
        df_supports  : pd.DataFrame,
        df_densities : pd.DataFrame, 
        lqd_sup=None,
        verbose_=True):
    """_summary_

    Args:
        df_supports (pd.DataFrame): support values for the densities
        df_densities (pd.DataFrame): densities corresponding to the supports
        lqd_sup_M (_type_, optional): support for the LQDT output. Defaults to [0,1] with nrows(df_densities)
            grid points.

    Returns:
        _type_: returns the LQDT (shared) support, the LQDT dataframe of values and the c=F(0) process

    Example: lqdSup, df_lqds, c = obtain_lqds(df_supports, df_densities) 
    """
    if lqd_sup is not None:
        M = lqd_sup.shape[0]
    else:
        M = df_densities.shape[0]
        lqd_sup = np.linspace(0,1, M)

    # lqds_sup = []
    lqds = []
    cs = []
    t0s = []
    cols = df_densities.columns
    for col in cols:
        density_support = df_supports.loc[:,col]
        density = df_densities.loc[:,col]
        t0_ = density_support.iloc[np.argmin(np.abs(density_support))] #closes value to 0
        lqdsup, lqd, c = dens2lqd(
                                dens = density, 
                                dSup = density_support, 
                                lqdSup=lqd_sup, 
                                t0 = t0_,
                                verbose=verbose_
                                )
        # lqds_sup.append(lqdSup) # no need since the image of the LQDT is shared by all densities by construction
        lqds.append(pd.Series(lqd))
        cs.append(c)
        t0s.append(t0_)
    
    df_lqds = pd.concat(lqds, axis=1)
    df_lqds.columns = cols
    
    return lqdsup, df_lqds, cs, t0s


def lqd2dens(
    lqd,
    lqdSup=None,
    dSup=None,
    t0=0.0,
    c=0.0,
    useSplines=True,
    cut=(0, 0),
    verbose=True
):
    """
    Reconstruct a probability density function from its Log–Quantile–Density
    (LQD) transform. This is the inverse of `dens2lqd`.

    This implementation follows the mathematical framework used in
    Kokoszka & Reimherr (2017, 2019) for functional transformations of 
    probability density functions.

    ----------------------------------------------------------------------
    Inverting the Transform
    ----------------------------------------------------------------------
    The goal is to reconstruct f from L.

    Step 1 — Recover the quantile density:
        q(u) = exp(L(u)).                                 (2)

    Step 2 — Recover the quantile function via integration:
        Q(u) = t0 + ∫_0^u q(s) ds  −  constant.            (3)

    The required constant is chosen so that
        Q(c) = t0,
    where c = F(t0).

    This ensures that the recovered quantile aligns with the required
    probability level c.

    Step 3 — The density is obtained by the identity:
        f(x) = 1 / q(F(x)).                                (4)

    Numerically, we compute dtemp = Q(u) and dens_temp = exp(-L(u)).
    Then we interpolate
        dens(x) = dens_temp(Q^{-1}(x)).

    Finally, density is renormalized to integrate to 1.

    Parameters
    ----------
    lqd : array-like
        Log quantile density on lqdSup.
    lqdSup : array-like, optional
        Support for lqd domain. Must start at 0, end at 1.
    dSup : unused (kept for compatibility)
    t0 : float
        Value such that the target CDF satisfies F(t0) = c.
    c : float
        CDF value in lqdSup at t0.
    useSplines : bool
        Whether to use splines when integrating exp(lqd).
    cut : tuple(int,int)
        How many boundary points to remove on left and right.
    verbose : bool
        Print messages?


    ----------------------------------------------------------------------
    Returns
    ----------------------------------------------------------------------
    dSup : ndarray
        Support grid for the reconstructed density. Equal to the 
        reconstructed quantile Q(u) after monotonicity correction.

    dens : ndarray
        Reconstructed density f(x). Renormalized to integrate to 1.
    """

    lqd = np.asarray(lqd)

    # --- handle lqdSup ---
    if lqdSup is None:
        lqdSup = np.linspace(0, 1, len(lqd))
    else:
        lqdSup = np.asarray(lqdSup)

    # --- Check support boundaries ---
    if not np.allclose([lqdSup.min(), lqdSup.max()], [0, 1]):
        if verbose:
            print("Problem with LQD domain boundaries — resetting to default.")
        lqdSup = np.linspace(0, 1, len(lqd))

    M = len(lqd)
    cut = list(cut)

    # --- find infinite exp(lqd) ---
    r = np.where(np.exp(lqd) == np.inf)[0]
    if len(r) > 0:
        mid = M // 2
        if np.any(r < mid):
            cut[0] = max(cut[0], r[r < mid].max())
        if np.any(r >= mid):
            cut[1] = max(cut[1], M - r[r >= mid].min() - 1)

    # --- cut boundaries ---
    lqdSup = lqdSup[cut[0]: M - cut[1]]
    lqd = lqd[cut[0]: M - cut[1]]
    M = len(lqd)

    # --- ensure c is in lqdSup ---
    if c not in lqdSup:
        if c < lqdSup[0] or c > lqdSup[-1]:
            raise ValueError("c is not within range of lqdSup after cutoff.")
        if verbose:
            print("c not in lqdSup — resetting to closest value.")
        c = lqdSup[np.argmin(np.abs(lqdSup - c))]

    c_ind = np.where(lqdSup == c)[0][0]

    # --- compute dtemp ---
    if useSplines:
        sp = InterpolatedUnivariateSpline(lqdSup, lqd, k=3)
        fexp = lambda t: np.exp(sp(t))

        integrals = [0]
        for i in range(1, M):
            v, _ = quad(fexp, lqdSup[i - 1], lqdSup[i])
            integrals.append(v)

        dtemp = t0 + np.cumsum(integrals)

        shift_v, _ = quad(fexp, lqdSup[0], lqdSup[c_ind])
        dtemp = dtemp - shift_v

    else:
        exp_lqd = np.exp(lqd)
        dtemp = t0 + cumulative_trapezoid(exp_lqd, lqdSup, initial=0)
        shift_v = np.trapezoid(exp_lqd[:c_ind + 1], lqdSup[:c_ind + 1])
        dtemp = dtemp - shift_v

    # =====================================================
    #     CORRECT R-LIKE DUPLICATE REMOVAL (NO ERRORS)
    # =====================================================

    mid = M // 2

    indL = duplicated_tol(dtemp[:mid], tol=1e-10, from_last=True)
    indR = duplicated_tol(dtemp[mid:], tol=1e-10, from_last=False)

    keep = ~(np.concatenate([indL, indR]))

    # ---- SAVE ORIGINAL ENDPOINTS (CRITICAL FIX) ----
    d0 = dtemp[0]
    d1 = dtemp[-1]

    dtemp = dtemp[keep]
    dens_temp = np.exp(-lqd[keep])

    # =====================================================
    # Interpolate & Normalize  (R faithful)
    # =====================================================

    # Build support from ORIGINAL endpoints, not shortened dtemp
    dSup = np.linspace(d0, d1, M)

    dens = np.interp(
        dSup,
        dtemp,
        dens_temp,
        left=dens_temp[0],
        right=dens_temp[-1]
    )

    # Normalize density
    area = np.trapezoid(dens, dSup)
    dens = dens / area * (lqdSup[-1] - lqdSup[0])
    
    return dSup, dens

def obtain_densities_from_lqd(
        df: pd.DataFrame, 
        lqdSup_: np.array, 
        c_: np.array,
        t_: np.array,
        cut_invalid: bool=True,
        verbose=True):
    """_summary_

    Args:
        df (pd.DataFrame): lqdensities dataframe
        lqdSup_ (np.array): common (shared) support of lqdensities
        c_ (np.array): F(0) for each lqdensity

    Returns:
        _type_: _description_
    """
    cols = df.columns

    supports  = []
    densities = []
    i=0
    for col in cols:
        lqd = df.loc[:, col]
        if cut_invalid:
            cut = compute_lqd_cut(lqd, cut_nan=True)
        else:
            cut = (0,0)
        backward_support, backward_density = lqd2dens(
                                                    lqd, 
                                                    lqdSup_, 
                                                    c = c_[i], 
                                                    t0=t_[i], 
                                                    cut=cut,
                                                    verbose=verbose
                                                    )
        supports.append(pd.Series(backward_support))
        densities.append(pd.Series(backward_density))
        i += 1

    df_backward_supports = pd.concat(supports, axis=1) 
    df_backward_supports.columns = cols
    df_backward_densities = pd.concat(densities, axis=1)
    df_backward_densities.columns = cols

    return df_backward_supports, df_backward_densities


class mLQDT:
    """
    bovespa_mLQDT = mLQDT(
                    df_densities,
                    df_densities_supports
                )
    bovespa_mLQDT.densities_to_lqdensities()
    """
    def __init__(
            self, 
            densities : pd.DataFrame,
            densities_supports : pd.DataFrame
            ):
        
        # Data
        self.Y         = densities
        self.Y_support = densities_supports

        # Metadata
        self.columns  = densities.columns
        self.Y_n_rows = densities.shape[0]
        self.Y_n_cols = densities.shape[1]

        # Results 
        self.lqd_support = None
        self.lqd         = None
        self.c           = None
        self.t           = None


    def densities_to_lqdensities(self, 
                                 lqd_support=None,
                                 verbose=True):
        lqdSup, df_lqds, c, t = obtain_lqds(
                                        self.Y_support, 
                                        self.Y,
                                        lqd_sup = lqd_support,
                                        verbose_=verbose)
        
        self.lqd_support = lqdSup
        self.lqd         = df_lqds
        self.c           = c
        self.t           = t

    def lqdensities_to_densities(self, 
                                 adhoc_lqd_support,
                                 adhoc_lqds,
                                 adhoc_c=None,
                                 verbose=True
                                ):
        
        if adhoc_c is not None:
            if verbose:
                print("'c' was not provided. Using original fitted values.")
            adhoc_c = self.c

        densities_supports, densities = obtain_densities_from_lqd(
                                                            adhoc_lqds,
                                                            adhoc_lqd_support,
                                                            c=adhoc_c,
                                                            verbose = True
                                        )
        
        return [densities_supports, densities]