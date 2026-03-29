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