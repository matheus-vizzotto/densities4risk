import numpy as np
import pandas as pd
from scipy.integrate import simpson, trapezoid
from scipy.interpolate import interp1d

def integrates_one(
        u_grid:  np.array,
        density: np.array,
        tol:     float = 1e-4
    ) -> bool:
    from scipy.integrate import simpson

    integral = simpson(y=density, x=u_grid)
    is_normalized = np.isclose(integral, 1.0, atol=tol)
    
    return is_normalized


def verify_density_area(x, f_x, method='simpson'):
    """
    Verifies the area under f(x) using the chosen numerical method.
    
    Args:
        x (array): The domain (e.g., Bitcoin returns).
        f_x (array): The density values.
        method (str): 'simpson' or 'trapezoid'.
    """
    if method.lower() == 'simpson':
        area = simpson(y=f_x, x=x)
    elif method.lower() == 'trapezoid':
        area = trapezoid(y=f_x, x=x)
    else:
        raise ValueError("Method must be 'simpson' or 'trapezoid'")
        
    print(f"[{method.capitalize()}] Total Area: {area:.8f}")
    return area

def _calculate_density_mean(
        u_grid: np.array, 
        density: np.array
        ) -> float:
    """
    Calculates the first raw moment: mu = integral(u * f(u) du)
    """
    mean = simpson(y=u_grid * density, x=u_grid)

    return mean

def _calculate_density_variance(
        u_grid: np.array, 
        density: np.array
        ) -> float:
    """
    Calculates the second central moment: Var = integral((u - mu)^2 * f(u) du)
    """
    mu = _calculate_density_mean(u_grid, density)
    
    # Calculate (u - mu)^2
    squared_deviation = (u_grid - mu)**2
    
    # Integrate squared deviation weighted by the density
    variance = simpson(y=squared_deviation * density, x=u_grid)
    
    return variance


def process_bootstrap_kdes(kde_list, grid_points=5001):
    """
    Interpolates a list of KDE DataFrames to a common grid, 
    returns a single aligned DataFrame and includes the pointwise mean.
    """
    if not kde_list:
        return None

    # 1. Establish the global support grid
    min_support = min(df.index.min() for df in kde_list)
    max_support = max(df.index.max() for df in kde_list)
    common_grid = np.linspace(min_support, max_support, num=grid_points)

    # 2. Interpolate each KDE and store raw values in a list
    # We use the original column names (e.g., '2026-05-12' from image_cb66c0.png)
    column_names = []
    all_values = []

    for i, df in enumerate(kde_list):
        df_sorted = df.sort_index()
        
        # Interpolate onto the common grid
        y_interp = np.interp(
            x=common_grid, 
            xp=df_sorted.index, 
            fp=df_sorted.iloc[:, 0], 
            left=0, 
            right=0
        )
        
        all_values.append(y_interp)
        # Use existing column name or fallback to a numbered bootstrap
        column_names.append(df.columns[0] if not df.columns.empty else f"boot_{i}")

    # 3. Create the master DataFrame
    # np.column_stack turns the list of arrays into a 2D matrix (rows x bootstraps)
    interpolated_df = pd.DataFrame(
        data=np.column_stack(all_values),
        index=common_grid,
        columns=column_names
    )

    # 4. Add the pointwise mean as a final column
    interpolated_df["pointwise_mean"] = interpolated_df.mean(axis=1)

    return interpolated_df


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