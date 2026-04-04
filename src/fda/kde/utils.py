import numpy as np
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