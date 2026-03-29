import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import pywt

from typing import List, Tuple

def inner_product(f, g, du):
    """
    Compute the L^2 inner product ⟨f, g⟩ = ∫ f(u) g(u) du
    in its discrete approximation.

    Parameters
    ----------
    f : ndarray (m, k1)
        First function(s) evaluated on a grid.
    g : ndarray (m, k2)
        Second function(s) evaluated on the same grid.
    du : float
        Grid spacing.

    Returns
    -------
    ndarray (k1, k2)
        Matrix of inner products.
    """
    return np.dot(f.T, g) * du


def L2norm(f, du):
    """
    Compute the L^2 norm of a function vector f.

    Parameters
    ----------
    f : ndarray (m,)
        Function values on a grid.
    du : float
        Grid spacing.

    Returns
    -------
    float
        The L^2 norm.
    """
    return np.sqrt(inner_product(f, f, du))



############################
######### Wavelets #########
############################

def wavedec_sizes(signal_length, wavelet_name, level):
    """
    Compute expected pywt.wavedec sizes:
      returns sizes = [len(A_N), len(D_N), ..., len(D_1)]
      and total = sum(sizes)
    """
    w = pywt.Wavelet(wavelet_name)
    F = w.dec_len
    L = signal_length

    lengths = []
    Lj = L
    for j in range(level):
        Lj = (Lj + F - 1) // 2
        lengths.append(Lj)

    sizes = [lengths[-1]] + lengths[::-1]

    total = sum(sizes)
    return sizes, total

# split wavelet vector into coeff lists
def vec2coeffs(vec, sizes):
    coeffs = []
    pos = 0
    for s in sizes:
        coeffs.append(vec[pos:pos+s])
        pos += s
    return coeffs


############################
######### GENERAL ##########
############################
def align_and_normalize_density(x_obs, f_obs, x_hat, f_hat,
                                x_common=None, n_points=5001,
                                interp_kind='linear', eps=1e-12):
    """
    Interpolate two density curves onto a common grid, clip negatives,
    and renormalize so each integrates to 1.

    Returns: x_common, f_obs_common, f_hat_common
    """

    # Choose a common grid covering both supports
    if x_common is None:
        lo = min(np.nanmin(x_obs), np.nanmin(x_hat))
        hi = max(np.nanmax(x_obs), np.nanmax(x_hat))
        x_common = np.linspace(lo, hi, n_points)

    # safe interpolation: fill outside range with 0
    interp_obs = interp1d(x_obs, f_obs, kind=interp_kind,
                          bounds_error=False, fill_value=0.0)
    interp_hat = interp1d(x_hat, f_hat, kind=interp_kind,
                          bounds_error=False, fill_value=0.0)

    f_obs_c = interp_obs(x_common)
    f_hat_c = interp_hat(x_common)

    # sanitize: remove NaN/inf and negatives
    f_obs_c = np.nan_to_num(f_obs_c, nan=0.0, posinf=0.0, neginf=0.0)
    f_hat_c = np.nan_to_num(f_hat_c, nan=0.0, posinf=0.0, neginf=0.0)
    f_obs_c = np.maximum(f_obs_c, 0.0)
    f_hat_c = np.maximum(f_hat_c, 0.0)

    # renormalize each to integrate to 1
    area_obs = np.trapezoid(f_obs_c, x_common)
    area_hat = np.trapezoid(f_hat_c, x_common)

    # avoid division by zero
    if area_obs <= eps:
        raise ValueError("Observed density integrates to ~0 after interpolation.")
    if area_hat <= eps:
        raise ValueError("Forecast density integrates to ~0 after interpolation.")

    f_obs_c = f_obs_c / area_obs
    f_hat_c = f_hat_c / area_hat

    return x_common.copy(), f_obs_c.copy(), f_hat_c.copy()


def align_densities(
        densities_support : pd.DataFrame,
        densities : pd.DataFrame,
        fc_densities_support : pd.DataFrame,
        fc_densities : pd.DataFrame,
        cols : List,
        n : int = 5001
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Align observed and forecast density functions onto a common grid
    at each time index and renormalize them for divergence evaluation.

    This function is designed for density forecast evaluation under
    time-varying supports, as in the modified log-quantile-density
    transformation (LQDT) framework of Kokoszka et al. (2019).
    For each selected column (typically corresponding to a forecast
    origin or time index), the observed and forecast densities are
    interpolated onto a common grid defined on the union of their
    supports. The resulting curves are truncated at zero and
    renormalized to integrate to one.

    Parameters
    ----------
    densities_support : pandas.DataFrame
        Support grids of the observed densities. Each column corresponds
        to a time index and contains the support on which the observed
        density is defined.

    densities : pandas.DataFrame
        Observed density values evaluated on `densities_support`.
        Must have the same shape and column labels as
        `densities_support`.

    fc_densities_support : pandas.DataFrame
        Support grids of the forecast densities. Each column corresponds
        to a time index and may differ from the observed support.

    fc_densities : pandas.DataFrame
        Forecast density values evaluated on `fc_densities_support`.
        Must have the same shape and column labels as
        `fc_densities_support`.

    cols : list
        List of column names (time indices) for which densities are
        aligned and normalized.

    Returns
    -------
    df_supp : pandas.DataFrame
        Common support grids obtained for each column after alignment.
        Each column contains the grid on which both observed and
        forecast densities are evaluated.

    df_f : pandas.DataFrame
        Observed densities interpolated onto the common grids and
        renormalized to integrate to one.

    df_fhat : pandas.DataFrame
        Forecast densities interpolated onto the common grids and
        renormalized to integrate to one.

    Notes
    -----
    * Alignment is performed independently for each column (time index),
      allowing for time-varying supports.
    * Interpolation, truncation at zero, and renormalization are handled
      by the auxiliary function `align_and_normalize_density`.
    * The returned densities are suitable for the computation of
      divergence-based accuracy measures such as the Kullback–Leibler
      divergence or Jensen–Shannon divergence.
    """

    base = pd.DataFrame(np.full([n,len(cols)], np.nan), columns = cols)

    df_supp = base.copy()
    df_f = base.copy()
    df_fhat = base.copy()
    for col in cols:
        supp_col, f_col, f_hat_col = align_and_normalize_density(densities_support.loc[:,col], densities.loc[:,col], fc_densities_support.loc[:,col], fc_densities.loc[:,col])
        df_supp.loc[:,col] = supp_col
        df_f.loc[:,col] = f_col
        df_fhat.loc[:,col] = f_hat_col
    
    return df_supp, df_f, df_fhat