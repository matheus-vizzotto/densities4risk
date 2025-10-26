# OBS: could not find function "trapzRcpp" -> install.packages("fdapace")

import numpy as np
import warnings
from scipy import integrate
from scipy import interpolate 


def T_q(
        dens: np.array, 
        dens_support: np.array, 
        lqd_support: np.array = None,
        t0 : int = None
        ) -> list[np.array, np.array, np.float64]:
    """Function for converting densities to log quantile density functions

    Args:
        dens (np.array): density values on dSup - must be strictly positive (otherwise will truncate) and integrate to 1.
        dens_support (np.array): support (grid) for Density domain.
        lqd_support (np.array, optional): support of length M for lqd domain - must begin at 0 and end at 1. Defaults to None.
        t0 (int, optional): value in dSup for which the cdf value c is retained, i.e. c = F(t0). Defaults to None.

    Returns:
    list[np.array, np.array, np.float64]: [lqd_support, lqd, c]
    """
    lqd_support = np.linspace(0,1, len(dens_support))
    t0 = dens_support[0]

    if np.any(dens<0):
        raise Exception("There are negative density values")
    if np.abs(np.trapezoid(dens, dens_support) - 1) > 1e-5:
        warnings.warn("Density does not integrate to 1 with tolerance of 1e-5 - renormalizing now.")
        dens = dens/np.trapezoid(dens, dens_support)
    if np.any(dens==0):  # FUNÇÃO NÃO CONSIDERA POSSIBILIDADE DE ZERO ENTRE LIMITES DE (SUPORTE>O)
        warnings.warn("There are some zero density values - truncating support grid so all are positive.")
        lower_bound = np.min(np.where(dens>0))
        upper_bound = np.max(np.where(dens>0))
        dens_support = dens_support[lower_bound : upper_bound]
        dens = dens/np.trapezoid(dens, dens_support) # ENTENDER POR QUE ISSO FUNCIONA
        
    N = len(dens_support)

    if lqd_support is None:
        lqd_support = np.linspace(0,1,N)
    elif not (np.array([lqd_support.min(),lqd_support.max()]) == np.array([0,1])).all():
        warnings.warn(("Problem with support of the LQD domain's boundaries - resetting to default."))
        lqd_support = np.linspace(0,1,N)
        
    if not (t0 in dens_support):
        warnings.warn("t0 is not a value in dSup - resetting to closest value")
        t0 = dens_support[np.argmin((abs(dens_support-t0)))]
    M = len(lqd_support)
    c_ind = np.where(dens_support == t0)[0][0]

    # Get CDF and lqd on temporary grid, compute c
    tmp = integrate.cumulative_trapezoid(dens, dens_support, initial=0)
    c = tmp[c_ind]

    N_tmp = len(tmp)

    # ---- Left duplicates (forward scan)
    first_half = tmp[:N_tmp//2]
    _, indL = np.unique(first_half, return_index=True)
    maskL = np.ones(len(first_half), dtype=bool)
    maskL[np.setdiff1d(np.arange(len(first_half)), indL)] = False  # mark duplicates as False

    # ---- Right duplicates (reverse scan)
    second_half = tmp[N_tmp//2:]
    _, indR_rev = np.unique(second_half[::-1], return_index=True)
    maskR = np.ones(len(second_half), dtype=bool)
    maskR[np.setdiff1d(np.arange(len(second_half)), len(second_half) - 1 - indR_rev)] = False

    mask = np.concatenate([maskL, maskR])
    qtemp = tmp[mask]

    lqd_temp = -np.log(dens[mask])

    # Interpolate lqd_support, keeping track of Inf values at boundary, then compute c
    lqd = np.full(M, 0)

    if np.any(np.isinf([lqd_temp[0], lqd_temp[N-1]])):
        tmpInd = np.arange(1, N)
        Ind = np.arange(1, M)
        if lqd_temp[0] == np.inf:
            lqd[0] = np.inf
            tmpInd = tmpInd[-1]
            Ind = Ind[-1]
        if lqd_temp[-1] == np.inf:
            lqd[M] = np.inf
            tmpInd = tmpInd[-len(tmpInd)]
            Ind = Ind[-len(Ind)]

        interp_values = np.interp(
            lqd_support[Ind-1],                 # xout
            qtemp[tmpInd-1],               # x
            lqd_temp[tmpInd-1],            # y
            left=lqd_temp[tmpInd-1][0],    # rule = 2 → use endpoint values
            right=lqd_temp[tmpInd-1][-1]
        )

        lqd[Ind] = interp_values
    else:
        interp_values = np.interp(
            lqd_support,
            qtemp,
            lqd_temp,
            left=lqd_temp[0],    # rule = 2 → use endpoint values
            right=lqd_temp[-1]
        )
        lqd = interp_values

    return [lqd_support, lqd, c]


def T_q_inv(
                lqd    : np.array, 
                lqdSup : np.array = None,
                t0 = 0,
                c  = 0,
                useSplines = True,
                cut = [0,0]
        ) -> list[np.array, np.array]:
        """Function for converting log quantile densities to densities.

        Args:
                lqd (np.array): log quantile density on lqdSup
                lqdSup (np.array, optional): support for lqd domain - must begin at 0 and end at 1. Defaults to None.
                t0 (int, optional): value for which the target cdf has F(t0) = c (default: 0). Defaults to 0.
                c (int, optional): value in lqdSup representing the value of target cdf at t0 (default: lqdSup[1]). Defaults to 0.
                useSplines (bool, optional): fit spline to the lqd when doing the numerical integration (default: TRUE). Defaults to True.
                cut (list, optional): vector with two elements, indicating how many boundary to cut on the left and right side (default: c(0, 0)).  More will be cut off if exp(lqd) is infinite for some values.. Defaults to [0,0].

        Returns:
                list[np.array, np.array]: [dSup, dens]
        """
        if not (np.array([lqdSup.min(),lqdSup.max()]) == np.array([0,1])).all():
                warnings.warn("Problem with support of the LQD domain's boundaries - resetting to default.")
                lqdSup = np.linspace(0,1,len(lqd))
        M = len(lqd)
        r = np.where(np.exp(lqd)==np.inf)[0]
        if len(r)>0:
                print(r)
                if np.any(r < int(M/2)):
                        cut[0] = np.max(cut[0], np.max(r[r<int(M/2)]))
                if np.any(r >= int(M/2)):
                        cut[1] = np.max(cut[1], M - np.min(r[r >= int(M/2)]))
        
        # Cut boundaries
        lqdSup = lqdSup[(cut[0]):(M-cut[1])]
        lqd    = lqd[(cut[0]):(M-cut[1])]
        M      = len(lqd) # reset N
        if not (c in lqdSup):
                if (c < lqdSup[0] or c > lqdSup[M]):
                        raise Exception ("c is not contained withing range of lqdSup after cutoff")
                print("c is not equal to a value in lqdSup - resetting to closest value")
                c = lqdSup[np.argmin(abs(lqdSup-c))]
        
        c_ind = np.where(lqdSup==c)[0]

        if useSplines:
                # 1. Natural cubic spline equivalent to splinefun(..., method="natural")
                lqd_sp = interpolate.CubicSpline(lqdSup, lqd, bc_type='natural')
                # 2. Define exp(lqd_sp(t))
                def lqd_exp(t):
                        return np.exp(lqd_sp(t))
                # 3. Compute the cumulative integrals between successive lqdSup points
                integrals = [integrate.quad(lqd_exp, lqdSup[i-1], lqdSup[i])[0] for i in range(1, len(lqdSup))]
                # 4. Construct dtemp (cumulative density grid)
                dtemp = t0 + np.concatenate(([0], np.cumsum(integrals))) - \
                                integrate.quad(lqd_exp, lqdSup[0], lqdSup[c_ind])[0]
        else:
                # Compute cumulative trapezoidal integral
                cum_int = integrate.cumulative_trapezoid(np.exp(lqd), lqdSup, initial=0)

                # Compute trapezoidal integral up to c_ind
                trap_int = integrate.trapezoid(np.exp(lqd[:c_ind]), lqdSup[:c_ind], initial=0)

                # Build dtemp
                dtemp = t0 + cum_int - trap_int
        # Remove Duplicates
        # --- Split indices
        mid = M // 2

        # --- Detect duplicates (left: check from end, right: from start)
        _, indL = np.unique(dtemp[:mid][::-1], return_index=True)
        indL_mask = np.ones(mid, dtype=bool)
        indL_mask[np.setdiff1d(np.arange(mid), mid - 1 - indL)] = False  # equivalent to duplicated(..., fromLast=TRUE)

        _, indR = np.unique(dtemp[mid:], return_index=True)
        indR_mask = np.ones(M - mid, dtype=bool)
        indR_mask[np.setdiff1d(np.arange(M - mid), indR)] = False  # equivalent to duplicated(...)

        # --- Combine masks
        mask = np.concatenate([indL_mask, indR_mask])

        # --- Apply mask
        dtemp = dtemp[mask]
        dens_temp = np.exp(-lqd[mask])

        # Interpolate to dSup and normalize
        # 1. Build uniform grid
        dSup = np.linspace(dtemp[0], dtemp[-1], M)

        # 2. Interpolate dens_temp to dSup (constant extrapolation on both sides)
        interp_func = interpolate.interp1d(dtemp, dens_temp, kind='linear',
                        bounds_error=False, fill_value=(dens_temp[0], dens_temp[-1]))
        dens = interp_func(dSup)

        # 3. Normalize, accounting for boundary cutoff
        area = integrate.trapezoid(dens, dSup)
        dens = dens / area * (lqdSup[-1] - lqdSup[0])

        return [dSup, dens]