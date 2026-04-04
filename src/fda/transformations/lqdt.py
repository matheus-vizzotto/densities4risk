import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from src.fda.transformations.utils import duplicated_tol_sorted, compute_lqd_cut 


# Fixed M to better approximate functions. Smaller M (like 256) generates worse approximations.
# M = 5001

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
    # Returns True if the value is either NaN or Inf (ADAPTATION FROM R)
    temp_first_inf = np.isnan(lqd_temp[0]) | np.isinf(lqd_temp[0])
    temp_last_inf = np.isnan(lqd_temp[-1]) | np.isinf(lqd_temp[-1])

    if temp_first_inf or temp_last_inf:
        
        print("NaNs detected.")

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

    # ---- cache exp(lqd) ----
    exp_lqd = np.exp(lqd)

    # --- find infinite exp(lqd) ---
    r = np.where(exp_lqd == np.inf)[0]
    if len(r) > 0:
        mid = M // 2
        if np.any(r < mid):
            cut[0] = max(cut[0], r[r < mid].max())
        if np.any(r >= mid):
            cut[1] = max(cut[1], M - r[r >= mid].min() - 1)

    # --- cut boundaries ---
    lqdSup = lqdSup[cut[0]: M - cut[1]]
    lqd = lqd[cut[0]: M - cut[1]]
    exp_lqd = exp_lqd[cut[0]: M - cut[1]]
    M = len(lqd)

    # --- ensure c is in lqdSup (closest index) ---
    c_ind = np.argmin(np.abs(lqdSup - c))
    c = lqdSup[c_ind]

    # --- compute dtemp ---
    if useSplines:
        sp = InterpolatedUnivariateSpline(lqdSup, lqd, k=3)
        q = np.exp(sp(lqdSup))

        Q = cumulative_trapezoid(q, lqdSup, initial=0)
        dtemp = t0 + Q - Q[c_ind]

    else:
        Q = cumulative_trapezoid(exp_lqd, lqdSup, initial=0)
        dtemp = t0 + Q - Q[c_ind]

    # =====================================================
    #     CORRECT R-LIKE DUPLICATE REMOVAL (NO ERRORS)
    # =====================================================

    mid = M // 2

    indL = duplicated_tol_sorted(dtemp[:mid], tol=1e-10, from_last=True)
    indR = duplicated_tol_sorted(dtemp[mid:], tol=1e-10, from_last=False)

    keep = ~(np.concatenate([indL, indR]))

    # ---- SAVE ORIGINAL ENDPOINTS (CRITICAL FIX) ----
    d0 = dtemp[0]
    d1 = dtemp[-1]

    dtemp = dtemp[keep]
    dens_temp = 1.0 / exp_lqd[keep]

    # =====================================================
    # Interpolate & Normalize  (R faithful)
    # =====================================================

    dSup = np.linspace(d0, d1, M)

    dens = np.interp(
        dSup,
        dtemp,
        dens_temp,
        left=dens_temp[0],
        right=dens_temp[-1]
    )

    area = np.trapezoid(dens, dSup)
    dens = dens / area * (lqdSup[-1] - lqdSup[0])

    return dSup, dens

#---------------------------------------------------------------------
#--------------------------- PIPELINES -------------------------------
#---------------------------------------------------------------------

def obtain_lqds(
        df_supports  : pd.DataFrame,
        df_densities : pd.DataFrame, 
        lqd_sup: np.array=None,
        fill_nan_absurd: float=20,
        verbose_=True
        ):
    """_summary_

    Args:
        df_supports (pd.DataFrame): support values for the densities
        df_densities (pd.DataFrame): densities corresponding to the supports
        lqd_sup_M (_type_, optional): support for the LQDT output. Defaults to [0,1] with nrows(df_densities)
            grid points.
        fill_nan_absurd (bool, optional): this argument substitutes NaN values at the boundaries with an extreme
            value that is to be cut by the lqd2dens function -- i.e., this only makes sense if lqd2dens outliers
            at the boundaries.

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
    t0 = 0
    cols = df_densities.columns
    for col in cols:
        density_support = df_supports.loc[:,col]
        density = df_densities.loc[:,col]
        t0 = density_support.iloc[np.argmin(np.abs(density_support))] #closes value to 0
        lqdsup, lqd, c = dens2lqd(
                                dens = density, 
                                dSup = density_support, 
                                lqdSup=lqd_sup, 
                                t0 = t0,
                                verbose=verbose_
                                )
        if fill_nan_absurd:
            if np.isnan(lqd[0]):
                lqd[0] = fill_nan_absurd
            if np.isnan(lqd[-1]):
                lqd[-1] = fill_nan_absurd

        # lqds_sup.append(lqdSup) # no need since the image of the LQDT is shared by all densities by construction
        lqds.append(pd.Series(lqd))
        cs.append(c)
        # t0s.append(t0_)
    
    df_lqds = pd.concat(lqds, axis=1)
    df_lqds.columns = cols
    
    return lqdsup, df_lqds, cs, t0

def obtain_densities_from_lqd(
        df: pd.DataFrame, 
        lqdSup_: np.array, 
        c_: np.array,
        t_: float,
        cut_invalid_boundaries: bool=True,
        ensure_integral_constraint: bool=True,
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
        if cut_invalid_boundaries:
            cut = compute_lqd_cut(lqd, cut_nan=True)
        else:
            cut = (0,0)
        backward_support, backward_density = lqd2dens(
                                                    lqd, 
                                                    lqdSup_, 
                                                    # c = c_[i], 
                                                    c=c_,
                                                    # t0=t_[i],
                                                    t0=t_, 
                                                    cut=cut,
                                                    verbose=verbose
                                                    )
        if ensure_integral_constraint:
            # print('ensuring integral constraint.')
            backward_density /= np.trapezoid(backward_density, x=backward_support, axis=0)
        supports.append(pd.Series(backward_support))
        densities.append(pd.Series(backward_density))
        i += 1

    df_backward_supports = pd.concat(supports, axis=1) 
    df_backward_supports.columns = cols
    df_backward_densities = pd.concat(densities, axis=1)
    df_backward_densities.columns = cols

    return df_backward_supports, df_backward_densities

class LQDRepresentation:
    """
    Container for the Log–Quantile–Density (LQD) representation of a
    collection of probability density functions.

    This object stores all curve-specific quantities required to
    reconstruct densities via the inverse LQD transform.

    ------------------------------------------------------------------
    Parameters
    ------------------------------------------------------------------
    lqd : pandas.DataFrame or ndarray
        LQD values for each curve evaluated on a common grid.
        Shape (M, T), where:
            M = number of grid points in [0,1]
            T = number of curves (e.g., time points)

    lqd_support : ndarray
        Grid in [0,1] where the LQD is evaluated.

    c : array-like
        Array of CDF values at the reference points t0 for each curve.
        Shape (T,).

    t0 : float
        Reference point in the original support such that
            c = F(t0)

    ------------------------------------------------------------------
    Attributes
    ------------------------------------------------------------------
    lqd : pandas.DataFrame or ndarray
    lqd_support : ndarray
    c : ndarray
    t0 : float

    ------------------------------------------------------------------
    Notes
    ------------------------------------------------------------------
    - This object is intentionally self-contained: all information needed
      for inverse transformation is stored here.
    """
    def __init__(self, lqd, lqd_support, c, t0):
        self.lqd         : pd.DataFrame = lqd
        self.lqd_support : np.array     = lqd_support
        self.c           : np.array     = c
        self.t0          : np.array     = t0
        self.du          : float        = lqd_support[1] - lqd_support[0]

class mLQDT:
    """
    Log–Quantile–Density Transformation (LQDT) operator.

    This class implements the forward and inverse transformations between
    probability density functions and their Log–Quantile–Density (LQD)
    representations, as described in Kokoszka (2019).

    ------------------------------------------------------------------
    Overview
    ------------------------------------------------------------------
    The transformation maps densities from a nonlinear space into an
    approximately linear Hilbert space:

        f(x)  →  L(u)

    enabling the application of Functional Data Analysis (FDA) tools.

    The inverse transformation reconstructs densities from their LQD
    representation.

    ------------------------------------------------------------------
    Methods
    ------------------------------------------------------------------
    transform(...)
        Converts a collection of densities into an LQDRepresentation.

    inverse_transform(...)
        Reconstructs densities from an LQDRepresentation.

    ------------------------------------------------------------------
    Design
    ------------------------------------------------------------------
    - All curve-specific quantities (LQD, c, t0) are stored in the
      LQDRepresentation object.
    - Supports DataFrame-based workflows where:
        columns = curves (e.g., dates)
        index   = support/grid values
    """
    def __init__(self):
        pass

    def transform(
            self,
            densities          : pd.DataFrame,
            densities_supports : pd.DataFrame,
            lqd_support        : np.array      = None,
            verbose            : bool          = True
            ):
        """
        Apply the Log–Quantile–Density transformation to a collection
        of densities.

        ------------------------------------------------------------------
        Parameters
        ------------------------------------------------------------------
        densities : pandas.DataFrame
            Density values. Each column corresponds to a curve.

        densities_supports : pandas.DataFrame
            Support values corresponding to each density curve.
            Must have the same shape and column structure as `densities`.

        lqd_support : ndarray, optional
            Grid in [0,1] where the LQD will be evaluated.
            If None, a default uniform grid is used.

        verbose : bool, default=True
            Whether to print warnings from the underlying transformation.

        ------------------------------------------------------------------
        Returns
        ------------------------------------------------------------------
        LQDRepresentation
            Object containing:
                - LQD values
                - LQD support grid
                - c values
                - reference point t0

        ------------------------------------------------------------------
        Notes
        ------------------------------------------------------------------
        - Each curve is transformed independently.
        - The output grid is shared across curves.
        """
        lqd_support, df_lqds, cs, t0 = obtain_lqds(
                                        densities_supports, 
                                        densities,
                                        lqd_sup = lqd_support,
                                        verbose_=verbose
                                        )

        return LQDRepresentation(
                            lqd         = df_lqds, 
                            lqd_support = lqd_support, 
                            c           = cs, 
                            t0          = t0
                            )

    def inverse_transform(
                        self, 
                        lqd_obj:        LQDRepresentation,
                        common_support: bool = False,
                        verbose:        bool = True
                        ):
        """
        Reconstruct probability density functions from their
        Log–Quantile–Density representation.

        ------------------------------------------------------------------
        Parameters
        ------------------------------------------------------------------
        lqd_obj : LQDRepresentation
            Object containing LQD values and all necessary metadata
            (lqd_support, c, t0).

        verbose : bool, default=True
            Whether to print warnings during reconstruction.

        ------------------------------------------------------------------
        Returns
        ------------------------------------------------------------------
        densities_supports : pandas.DataFrame or ndarray
            Reconstructed support grids for each curve.

        densities : pandas.DataFrame or ndarray
            Reconstructed density values.

        ------------------------------------------------------------------
        Notes
        ------------------------------------------------------------------
        - Reconstruction is performed independently for each curve.
        - The resulting densities are normalized to integrate to 1.
        - The original support is recovered via the quantile function.
        """
        lqds = lqd_obj.lqd
        lqd_support = lqd_obj.lqd_support
        cs = lqd_obj.c
        t0s = lqd_obj.t0

        densities_supports, densities = obtain_densities_from_lqd(
                                                            df      =   lqds,
                                                            lqdSup_ =   lqd_support,
                                                            c_      =   cs,
                                                            t_      =   t0s,
                                                            verbose = verbose
                                        )

        if not common_support:
            return densities_supports, densities

        # 2. Define the global support range
        grid_size = len(lqd_support)
        global_min = densities_supports.min().min()
        global_max = densities_supports.max().max()
        interpolated_support = np.linspace(global_min, global_max, grid_size)

        # 3. Interpolate each curve
        interp_df = pd.DataFrame(index=interpolated_support, columns=densities.columns)
        
        for col in densities.columns:
            # Linear interpolation with zero-padding for tails
            f_interp = interp1d(
                densities_supports[col], 
                densities[col], 
                kind='linear', 
                bounds_error=False, 
                fill_value=0
            )
            
            dens_values = f_interp(interpolated_support)
            
            # Re-normalize to ensure the integral is 1 on the new grid
            area = np.trapezoid(dens_values, interpolated_support)
            if area > 0:
                dens_values /= area
            
            interp_df[col] = dens_values
            
        return interpolated_support, interp_df