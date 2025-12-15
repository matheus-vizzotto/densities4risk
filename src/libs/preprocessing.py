import numpy as np
from scipy.interpolate import interp1d, PchipInterpolator
import matplotlib.pyplot as plt

def align_and_normalize_density(x_obs, f_obs, x_hat, f_hat,
                                x_common=None, n_points=1000,
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
    area_obs = np.trapz(f_obs_c, x_common)
    area_hat = np.trapz(f_hat_c, x_common)

    # avoid division by zero
    if area_obs <= eps:
        raise ValueError("Observed density integrates to ~0 after interpolation.")
    if area_hat <= eps:
        raise ValueError("Forecast density integrates to ~0 after interpolation.")

    f_obs_c = f_obs_c / area_obs
    f_hat_c = f_hat_c / area_hat

    return x_common, f_obs_c, f_hat_c