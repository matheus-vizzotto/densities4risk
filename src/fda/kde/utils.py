import numpy as np
import pandas as pd
from scipy.integrate import simpson, trapezoid

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