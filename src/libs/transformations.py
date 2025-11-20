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


# def T_q(
#         density: np.array, 
#         density_support: np.array, 
#         lqd_support: np.array = None,
#         t0 : int = None
#         ) -> list[np.array, np.array, np.float64]:
#     """Function for converting densities to log quantile density functions

#     Args:
#         density (np.array): density values on density_support - must be strictly positive (otherwise will truncate) and integrate to 1.
#         density_support (np.array): support (grid) for Density domain.
#         lqd_support (np.array, optional): support of length M for lqd domain - must begin at 0 and end at 1. Defaults to None.
#         t0 (int, optional): value in density_support for which the cdf value c is retained, i.e. c = F(t0). Defaults to None.

#     Returns:
#         list[np.array, np.array, np.float64]: [lqd_support, lqd, c]

#     Example:
#         dens_support = np.linspace(0, 1, 512)
#         dens = np.full(512, 0.5)
#         lqd_support, lqd, c = T_q(density=dens, density_support=dens_support)
#     """
#     lqd_support = np.linspace(0,1, len(density_support))
#     t0 = density_support[0]

#     if np.any(density<0):
#         raise Exception("There are negative density values")
#     if np.abs(np.trapezoid(density, density_support) - 1) > 1e-5:
#         warnings.warn("Density does not integrate to 1 with tolerance of 1e-5 - renormalizing now.")
#         density = density/np.trapezoid(density, density_support)
#     if np.any(density==0):  
#         warnings.warn("There are some zero density values - truncating support grid so all are positive.")
#         lower_bound = np.min(np.where(density>0))
#         upper_bound = np.max(np.where(density>0))
#         density_support = density_support[lower_bound : upper_bound]
#         density = density/np.trapezoid(density, density_support) 
        
#     N = len(density_support)

#     if lqd_support is None:
#         lqd_support = np.linspace(0,1,N)
#     elif not (np.array([lqd_support.min(),lqd_support.max()]) == np.array([0,1])).all():
#         warnings.warn(("Problem with support of the LQD domain's boundaries - resetting to default."))
#         lqd_support = np.linspace(0,1,N)
        
#     if not (t0 in density_support):
#         warnings.warn("t0 is not a value in density_support - resetting to closest value")
#         t0 = density_support[np.argmin((abs(density_support-t0)))]
#     M = len(lqd_support)
#     c_ind = np.where(density_support == t0)[0][0]

#     # Get CDF and lqd on temporary grid, compute c
#     tmp = integrate.cumulative_trapezoid(density, density_support, initial=0)
#     c = tmp[c_ind]

#     N_tmp = len(tmp)

#     # ---- Left duplicates (forward scan)
#     first_half = tmp[:N_tmp//2]
#     _, indL = np.unique(first_half, return_index=True)
#     maskL = np.ones(len(first_half), dtype=bool)
#     maskL[np.setdiff1d(np.arange(len(first_half)), indL)] = False  # mark duplicates as False

#     # ---- Right duplicates (reverse scan)
#     second_half = tmp[N_tmp//2:]
#     _, indR_rev = np.unique(second_half[::-1], return_index=True)
#     maskR = np.ones(len(second_half), dtype=bool)
#     maskR[np.setdiff1d(np.arange(len(second_half)), len(second_half) - 1 - indR_rev)] = False

#     mask = np.concatenate([maskL, maskR])
#     qtemp = tmp[mask]

#     lqd_temp = -np.log(density[mask])

#     # Interpolate lqd_support, keeping track of Inf values at boundary, then compute c
#     lqd = np.full(M, 0)

#     if np.any(np.isinf([lqd_temp[0], lqd_temp[N-1]])):
#         tmpInd = np.arange(1, N)
#         Ind = np.arange(1, M)
#         if lqd_temp[0] == np.inf:
#             lqd[0] = np.inf
#             tmpInd = tmpInd[-1]
#             Ind = Ind[-1]
#         if lqd_temp[-1] == np.inf:
#             lqd[M] = np.inf
#             tmpInd = tmpInd[-len(tmpInd)]
#             Ind = Ind[-len(Ind)]

#         interp_values = np.interp(
#             lqd_support[Ind-1],            # xout
#             qtemp[tmpInd-1],               # x
#             lqd_temp[tmpInd-1],            # y
#             left=lqd_temp[tmpInd-1][0],    # rule = 2 → use endpoint values
#             right=lqd_temp[tmpInd-1][-1]
#         )

#         lqd[Ind] = interp_values
#     else:
#         interp_values = np.interp(
#             lqd_support,
#             qtemp,
#             lqd_temp,
#             left=lqd_temp[0],    # rule = 2 → use endpoint values
#             right=lqd_temp[-1]
#         )
#         lqd = interp_values

#     return [lqd_support, lqd, c]


# def T_q_inv(
#                 lqd    : np.array, 
#                 lqd_support : np.array = None,
#                 t0 = 0,
#                 c  = 0,
#                 useSplines = True,
#                 cut = [0,0]
#         ) -> list[np.array, np.array]:
#         """Function for converting log quantile densities to densities.

#         Args:
#             lqd (np.array): log quantile density on lqd_support
#             lqd_support (np.array, optional): support for lqd domain - must begin at 0 and end at 1. Defaults to None.
#             t0 (int, optional): value for which the target cdf has F(t0) = c (default: 0). Defaults to 0.
#             c (int, optional): value in lqd_support representing the value of target cdf at t0 (default: lqd_support[1]). Defaults to 0.
#             useSplines (bool, optional): fit spline to the lqd when doing the numerical integration (default: TRUE). Defaults to True.
#             cut (list, optional): vector with two elements, indicating how many boundary to cut on the left and right side (default: c(0, 0)).  More will be cut off if exp(lqd) is infinite for some values.. Defaults to [0,0].

#         Returns:
#             list[np.array, np.array]: [density_support, density]

#         Example:
#             y_lqd = np.full(512, np.log(2))
#             lqd_sup = np.linspace(0,1,len(y_lqd))
#             density_support, density = T_q_inv(lqd=y_lqd, lqd_support=lqd_sup)
#         """
#         if not (np.array([lqd_support.min(),lqd_support.max()]) == np.array([0,1])).all():
#                 warnings.warn("Problem with support of the LQD domain's boundaries - resetting to default.")
#                 lqd_support = np.linspace(0,1,len(lqd))
#         M = len(lqd)
#         r = np.where(np.exp(lqd)==np.inf)[0]
#         if len(r)>0:
#                 print(r)
#                 if np.any(r < int(M/2)):
#                         cut[0] = np.max(cut[0], np.max(r[r<int(M/2)]))
#                 if np.any(r >= int(M/2)):
#                         cut[1] = np.max(cut[1], M - np.min(r[r >= int(M/2)]))
        
#         # Cut boundaries
#         lqd_support = lqd_support[(cut[0]):(M-cut[1])]
#         lqd    = lqd[(cut[0]):(M-cut[1])]
#         M      = len(lqd) # reset N
#         if not (c in lqd_support):
#                 if (c < lqd_support[0] or c > lqd_support[M]):
#                         raise Exception ("c is not contained withing range of lqd_support after cutoff")
#                 print("c is not equal to a value in lqd_support - resetting to closest value")
#                 c = lqd_support[np.argmin(abs(lqd_support-c))]
        
#         c_ind = np.where(lqd_support==c)[0]

#         if useSplines:
#                 # 1. Natural cubic spline equivalent to splinefun(..., method="natural")
#                 lqd_sp = interpolate.CubicSpline(lqd_support, lqd, bc_type='natural')
#                 # 2. Define exp(lqd_sp(t))
#                 def lqd_exp(t):
#                         return np.exp(lqd_sp(t))
#                 # 3. Compute the cumulative integrals between successive lqd_support points
#                 integrals = [integrate.quad(lqd_exp, lqd_support[i-1], lqd_support[i])[0] for i in range(1, len(lqd_support))]
#                 # 4. Construct dtemp (cumulative density grid)
#                 dtemp = t0 + np.concatenate(([0], np.cumsum(integrals))) - \
#                                 integrate.quad(lqd_exp, lqd_support[0], lqd_support[c_ind])[0]
#         else:
#                 # Compute cumulative trapezoidal integral
#                 cum_int = integrate.cumulative_trapezoid(np.exp(lqd), lqd_support, initial=0)

#                 # Compute trapezoidal integral up to c_ind
#                 trap_int = integrate.trapezoid(np.exp(lqd[:c_ind]), lqd_support[:c_ind], initial=0)

#                 # Build dtemp
#                 dtemp = t0 + cum_int - trap_int
#         # Remove Duplicates
#         # --- Split indices
#         mid = M // 2

#         # --- Detect duplicates (left: check from end, right: from start)
#         _, indL = np.unique(dtemp[:mid][::-1], return_index=True)
#         indL_mask = np.ones(mid, dtype=bool)
#         indL_mask[np.setdiff1d(np.arange(mid), mid - 1 - indL)] = False  # equivalent to duplicated(..., fromLast=TRUE)

#         _, indR = np.unique(dtemp[mid:], return_index=True)
#         indR_mask = np.ones(M - mid, dtype=bool)
#         indR_mask[np.setdiff1d(np.arange(M - mid), indR)] = False  # equivalent to duplicated(...)

#         # --- Combine masks
#         mask = np.concatenate([indL_mask, indR_mask])

#         # --- Apply mask
#         dtemp = dtemp[mask]
#         dens_temp = np.exp(-lqd[mask])

#         # Interpolate to density_support and normalize
#         # 1. Build uniform grid
#         density_support = np.linspace(dtemp[0], dtemp[-1], M)

#         # 2. Interpolate dens_temp to density_support (constant extrapolation on both sides)
#         interp_func = interpolate.interp1d(dtemp, dens_temp, kind='linear',
#                         bounds_error=False, fill_value=(dens_temp[0], dens_temp[-1]))
#         density = interp_func(density_support)

#         # 3. Normalize, accounting for boundary cutoff
#         area = integrate.trapezoid(density, density_support)
#         density = density / area * (lqd_support[-1] - lqd_support[0])

#         return [density_support, density]





# Raw data to densities
def obtain_densities(
        df: pd.DataFrame, 
        M: int, 
        common_density_support=False):
    """Transform an nxT matrix (dataframe) into an MxT densities matrix, where
        n: number of observations in each input column (period)
        T: number of functional data objects
        M: output number of grid points where densities are estimated

    Args:
        df (pd.DataFrame): nxT matrix
        M (int): grid points for density evaluation
        common_density_support (bool, optional): Uses a common grid to evaluate densities by taking
        global minimum and maximum -- this is not advised for the LQDT since it can generate infinite
        boundary values because of values very close to zero. Defaults to False.

    Returns:
        _type_: returns one dataframe of supports and another of estimated densities, where the order
        of the columns on each are synchronized.

    Example: df_supports, df_densities = obtain_densities(df_returns, M=3000)
    """
    cols = df.columns

    if common_density_support:
        # 1) Global support
        global_min = df.min().min()
        global_max = df.max().max()
        u = np.linspace(global_min, global_max, M)

        # 2) Prepare density matrix (m × T)
        df_densities = pd.DataFrame(index=u, columns=cols)

        # 3) KDE for each day evaluated on a common support
        for t in cols:
            kde = gaussian_kde(df[t])
            df_densities[t] = kde(u)

        df_supports = u

    if not common_density_support:
        supports  = []
        densities = []
        for col in cols:
            data = df.loc[:, col]
            kde = gaussian_kde(data)
            left_endpoint = data.min()
            right_endpoint = data.max()
            x_grid = np.linspace(left_endpoint, right_endpoint, M)
            supports.append(pd.Series(x_grid))
            y_kde = kde(x_grid)
            densities.append(pd.Series(y_kde))
        df_supports = pd.concat(supports, axis=1)
        df_supports.columns = cols
        df_densities = pd.concat(densities, axis=1)
        df_densities.columns = cols

    return df_supports, df_densities



def dens2lqd(dens, dSup, lqdSup=None, t0=None, verbose=True):
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
        lqd_sup_M=None):
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
    if lqd_sup_M:
        M = lqd_sup_M
    else:
        M = df_densities.shape[0]
    lqdSup_ = np.linspace(0,1, M)

    # lqds_sup = []
    lqds = []
    cs = []
    cols = df_densities.columns
    for col in cols:
        density_support = df_supports.loc[:,col]
        density = df_densities.loc[:,col]
        t0_ = density_support.iloc[np.argmin(np.abs(density_support))] #closes value to 0
        lqdsup, lqd, c = dens2lqd(
                                dens = density, 
                                dSup = density_support, 
                                lqdSup=lqdSup_, 
                                t0 = t0_)
        # lqds_sup.append(lqdSup) # no need since the image of the LQDT is shared by all densities by construction
        lqds.append(pd.Series(lqd))
        cs.append(c)
    
    df_lqds = pd.concat(lqds, axis=1)
    df_lqds.columns = cols
    
    return lqdSup_, df_lqds, cs


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
    Python translation of the R function `lqd2dens`.

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

    Returns
    -------
    dict with:
        - dSup : density support
        - dens : density values
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

    # ----- Left half: duplicated(..., fromLast = TRUE) -----
    left = dtemp[:mid]
    indL = np.zeros(len(left), dtype=bool)
    seen = set()
    for i in range(len(left) - 1, -1, -1):
        if left[i] in seen:
            indL[i] = True
        seen.add(left[i])

    # ----- Right half: duplicated(..., fromLast = FALSE) -----
    right = dtemp[mid:]
    indR = np.zeros(len(right), dtype=bool)
    seen = set()
    for i in range(len(right)):
        if right[i] in seen:
            indR[i] = True
        seen.add(right[i])

    # Combine to full-length M mask
    keep = ~(np.concatenate([indL, indR]))
    dtemp = dtemp[keep]
    dens_temp = np.exp(-lqd[keep])

    # =====================================================
    #              Interpolate & Normalize
    # =====================================================

    dSup = np.linspace(dtemp[0], dtemp[-1], len(dtemp))
    dens = np.interp(dSup, dtemp, dens_temp)

    # Normalize density to integrate to 1 * boundary length
    area = np.trapezoid(dens, dSup)
    dens = dens / area * (lqdSup[-1] - lqdSup[0])

    return dSup, dens

def obtain_densities_from_lqd(
        df : pd.DataFrame, 
        lqdSup_ : np.array, 
        c_ : np.array):
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
        backward_support, backward_density = lqd2dens(lqd, lqdSup_, c = c_[i])
        supports.append(pd.Series(backward_support))
        densities.append(pd.Series(backward_density))
        i += 1

    df_backward_supports = pd.concat(supports, axis=1) 
    df_backward_densities = pd.concat(densities, axis=1)

    return df_backward_supports, df_backward_densities