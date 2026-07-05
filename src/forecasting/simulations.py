import numpy                        as np
import pandas                       as pd
import scipy.stats                  as stats
from scipy.stats import t
from scipy.special import gamma
import plotly.express               as px
import src.fda.transformations.lqdt as lqdt
import src.fda.plots                as fplt
import src.fda.simulations          as sim

def generate_ar_coeffs(d):
    """
    Generate AR coefficients for eigenfunctions dependence.
    """
    c = []
    for l0 in range(d):
        l = l0 + 1
        c0 = (-1)**l * (0.9 - ((0.5*l)/d))
        # print(c0, "\n")
        c.append(c0)
    return c

# def generate_ar_coeffs(K, p=1):
#     """
#     Generate stable AR(p) coefficients with decay across lags and components.
#     """
#     coeffs = np.zeros((K, p))

#     for k in range(K):
#         for j in range(p):
#             l = k + 1
#             coeffs[k, j] = ((-1)**(l + j)) * (0.8 / (j + 1)) * (1 - l/(2*K))

#     return coeffs

def generate_base_density(grid, kind="gaussian", **params):
    """
    Generates a normalized density f(x) evaluated on a grid using 
    flexible distribution parameters.

    Parameters
    ----------
    grid : ndarray
        The support points where the density is evaluated.
    kind : str
        The distribution type: 'gaussian', 'student_t', or 'mixture'.
    **params : dict
        Arbitrary keyword arguments passed to the scipy.stats distribution.
        Examples:
            - gaussian: mu=0, sigma=1 (mapped to loc/scale)
            - student_t: df=5, loc=0, scale=1
            - mixture: list of dicts under 'components'

    Examples
    --------
    >>> # 1. Standard Gaussian
    >>> f_norm = generate_base_density(grid, kind="gaussian", mu=0, sigma=1)
    
    >>> # 2. Heavy-Tailed Student-t 
    >>> f_heavy = generate_base_density(grid, kind="student_t", df=2.5, loc=0, scale=1)
    
    >>> # 4. Bimodal Mixture 
    >>> mix_params = {
    ...     "components": [
    ...         {"type": "norm", "weight": 0.4, "loc": -2, "scale": 0.8},
    ...         {"type": "norm", "weight": 0.6, "loc": 2, "scale": 1.2}
    ...     ]
    ... }
    >>> f_mix = generate_base_density(grid, kind="mixture", **mix_params)
    """
    if kind == "gaussian":
        # Map user-friendly names to scipy standard loc/scale
        loc = params.pop("mu", params.get("loc", 0))
        scale = params.pop("sigma", params.get("scale", 1))
        f = stats.norm.pdf(grid, loc=loc, scale=scale, **params)
        
    elif kind == "student_t":
        # Allows df, loc, scale, and even 'moments' or other scipy args
        f = stats.t.pdf(grid, **params)
        
    elif kind == "mixture":
        # Expects a list of dicts: [{'type': 'norm', 'weight': 0.5, 'mu': 0...}, ...]
        components = params.get("components", [
            {"type": "norm", "weight": 0.5, "loc": -2, "scale": 1},
            {"type": "norm", "weight": 0.5, "loc": 2, "scale": 1}
        ])
        
        f = np.zeros_like(grid)
        for comp in components:
            w = comp.pop("weight", 1/len(components))
            dist_type = comp.pop("type", "norm")
            # Dynamically call the distribution from scipy.stats
            dist_func = getattr(stats, dist_type)
            f += w * dist_func.pdf(grid, **comp)
            
    else:
        # Fallback: try to find the distribution directly in scipy.stats
        try:
            dist_func = getattr(stats, kind)
            f = dist_func.pdf(grid, **params)
        except AttributeError:
            raise ValueError(f"Distribution '{kind}' not recognized in scipy.stats")
    
    # Final normalization check
    total_area = np.trapezoid(f, grid)
    return f / total_area if total_area > 0 else f

def transform_to_l2(f, grid):
    """
    Apply modified Log Quantile Density Transform (mLQDT).
    """
    # Pseudo-steps (adapt to your implementation)
    F = np.cumsum(f) / np.sum(f)
    Q = np.interp(grid, F, grid)  # quantile function
    
    lqd = np.log(np.gradient(Q, grid))
    
    return lqd

def simulate_ar_scores(n, phis, sigma):
    """
    Simulates AR(1) score processes W_{kt}.
    
    Returns:
        scores: (K, n)
    """
    K = len(phis)
    scores = np.zeros((K, n))
    
    for k in range(K):
        shocks = np.random.normal(scale=sigma, size=n)
        for t in range(1, n):
            scores[k, t] = phis[k] * scores[k, t-1] + shocks[t]
    
    return scores

def build_basis(K, u, basis_type="cosine"):
    """
    Constructs an orthonormal trigonometric basis in the L2[0,1] Hilbert space.

    These basis functions are used to represent the stochastic components of 
    the Functional Time Series (FTS), such as the latent process X_t(u) or 
    the error term epsilon_t(u).

    Parameters
    ----------
    K : int
        Number of basis functions to generate (the dimension of the subspace).
    u : ndarray
        Discretization grid in the unit interval [0, 1]. Shape (M,).
    basis_type : {'cosine', 'sine'}, default='cosine'
        The type of trigonometric functions to use. 
        - 'cosine': {sqrt(2) * cos(j * pi * u)} for j = 1, ..., K
        - 'sine': {sqrt(2) * sin(j * pi * u)} for j = 1, ..., K

    Returns
    -------
    basis : ndarray
        An array of shape (K, M) containing the evaluated basis functions.

    Notes
    -----
    The basis is normalized such that the L2 norm of each function is 1, 
    assuming a uniform grid on [0,1]. Specifically:
    <phi_j, phi_l> = integral_{0}^{1} phi_j(u)phi_l(u)du = delta_{jl}
    """
    u = np.asarray(u)
    M = len(u)
    basis = np.zeros((K, M))

    for k in range(K):
        j = k + 1  # Standard indexing starting from 1
        if basis_type == "cosine":
            # Orthonormal Cosine Basis
            basis[k, :] = np.sqrt(2) * np.cos(j * np.pi * u)
        elif basis_type == "sine":
            # Orthonormal Sine Basis
            basis[k, :] = np.sqrt(2) * np.sin(j * np.pi * u)
        else:
            raise ValueError("basis_type must be 'sine' or 'cosine'")

    return basis

def build_dependence(scores, basis):
    """
    Computes X_t(u) = sum_k W_{kt} * phi_k(u)
    
    Returns:
        X: (nt, n)
    """
    return (scores.T @ basis).T

def simulate_brownian_bridge(n, nt, u, sigma_bb=0.05):
    """
    Simulates Brownian bridge noise ε_t(u).
    
    Returns:
        BB: (nt, n)
    """
    BB = np.zeros((nt, n))
    dt = 1 / (nt - 1)

    for t in range(n):
        dw = np.random.normal(scale=np.sqrt(dt), size=nt)
        w = np.cumsum(dw)
        bb = w - u * w[-1]
        BB[:, t] = sigma_bb * bb   # ← key change

    return BB

def simulate_curve_noise(n, nt, u, scale=1):
    """
    Simula ruído funcional usando uma combinação linear de senos.
    Baseado em Bhatia et al (2010).
    """
    mEps = np.zeros((nt, n))
    for ii in range(n):
        for jj in range(1, 11): 
            mEps[:, ii] += (
                stats.norm.rvs(scale=scale) * np.sqrt(2) *
                np.cos(np.pi * u * jj) / (2 ** (jj - 1))
            )

    return mEps

def simulate_white_noise(n, nt, sigma):
    """
    Simple iid functional noise.
    """
    return np.random.normal(scale=sigma, size=(nt, n))

def assemble_l2_process(base_lqd, X, noise):
    """
    Combines all components:
    
    Y_t(u) = base_lqd + X_t(u) + ε_t(u)
    """
    if base_lqd.ndim == 1:
        base_lqd = base_lqd[:, None]
        
    return base_lqd + X + noise

def simulate_l2_process(
    n,
    base_lqd,
    u,
    phis=None,
    K=None,
    sigma=0.05,
    error_sigma=0.05,
    noise_type="bridge",
    basis="cosine",
    as_dataframe=False,
    date_init=None,
    date_end=None
):
    """
    Simulate a functional time series in L2 via basis-driven dynamics.

    Generates curves Y_t(u) = base_lqd(u) + X_t(u) + ε_t(u), where X_t(u)
    is constructed from AR(1) score processes and an orthonormal cosine basis,
    and ε_t(u) is either Brownian bridge or white noise.

    Parameters
    ----------
    n : int
        Number of curves (time points).
    base_lqd : array-like, shape (nt,)
        Mean function in L2 space.
    u : array-like, shape (nt,)
        Grid over [0, 1].
    phis : array-like, shape (K,)
        AR(1) coefficients.
    sigma : float, optional
        Innovation standard deviation.
    noise_type : {"bb", "functional_noise", "white"}, optional
        Type of functional noise.
    as_dataframe : bool, optional
        If True, return Y as a pandas DataFrame with a daily date index.
    date_init : str or datetime-like, optional
        Start date for the index.
    date_end : str or datetime-like, optional
        End date for the index.

    Returns
    -------
    Y : ndarray or DataFrame, shape (nt, n)
        Simulated functional data.
    scores : ndarray, shape (K, n)
    X : ndarray, shape (nt, n)
    noise : ndarray, shape (nt, n)

    Notes
    -----
    If `as_dataframe=True` and no dates are provided, columns are set to the
    last `n` calendar days ending today.

    Example
    -------
    >>> Y, *_ = simulate_l2_process(
    ...     n=30,
    ...     base_lqd=base_lqd,
    ...     u=u,
    ...     phis=[0.8, 0.5],
    ...     as_dataframe=True
    ... )
    """

    nt = len(u)
    if phis is None:
        if K is None:
            raise ValueError("Provide either 'phis' or 'K'")
        phis = generate_ar_coeffs(K)
    else:
        phis = np.asarray(phis)
        K = len(phis)

    scores = simulate_ar_scores(n, phis, sigma)
    basis = build_basis(K, u, basis_type=basis)
    X = build_dependence(scores, basis)

    if noise_type == "bb":
        noise = simulate_brownian_bridge(n, nt, u, sigma_bb=sigma)
    elif noise_type == "functional_noise":
        noise = simulate_curve_noise(n, nt, u, error_sigma)
    elif noise_type == "white":
        noise = simulate_white_noise(n, nt, sigma)
    elif noise_type == "Null":
        noise = np.zeros_like(X)
    else:
        raise ValueError("noise_type must be 'bb', 'functional_noise' or 'white'")

    Y = assemble_l2_process(base_lqd, X, noise)

    if as_dataframe:
        # --- Date logic ---
        if date_init is None and date_end is None:
            dates = pd.date_range(end=pd.Timestamp.today().normalize(), periods=n, freq="D")
        elif date_end is not None and date_init is None:
            dates = pd.date_range(end=date_end, periods=n, freq="D")
        elif date_init is not None and date_end is None:
            dates = pd.date_range(start=date_init, periods=n, freq="D")
        else:
            dates = pd.date_range(start=date_init, end=date_end, freq="D")
            if len(dates) != n:
                raise ValueError("date range length must match n")

        Y = pd.DataFrame(Y, index=u, columns=dates)

    return Y, scores, X, noise, phis

class FDFSimulator:
    """
    Simulation engine for Functional Density Forecasting (FDF).
    
    Coordinates the transformation, simulation, and reconstruction of 
    density-valued time series.
    """
    def __init__(self, x_grid, u_grid=None, t0=None):
        self.x = x_grid  # Density support (e.g., returns)
        self.u = u_grid if u_grid is not None else np.linspace(0, 1, len(x_grid)) # LQD support [0, 1]
        # Default t0 to the point closest to zero
        self.t0 = t0 if t0 is not None else x_grid[np.argmin(np.abs(x_grid))]
        self.mlqdt = lqdt.mLQDT() 

        # fitted values
        self.result = None
        self.ar_coeffs = None

    def run_simulation(self, 
                       base_pdf, 
                       sigma=0.05, 
                       error_sigma=0.05,
                       dimensions=2,
                       phis=None, 
                       noise_type="bridge",
                       basis="cosine", 
                       n_curves=100, 
                       common_support=True
                       ):
        """
        Runs a full simulation pipeline using a user-provided baseline density.

        Parameters
        ----------
        base_pdf : np.ndarray
            The density values f(x) evaluated on self.x_grid.
        n_curves : int
            Number of functional observations to simulate.
        phis : list
            Autoregressive coefficients for the latent process.
        K: int
            Number of dimensions for the simulated process.
        sigma : float
            Scale of the functional noise.
        noise_type : str
            Type of stochastic perturbation ('bridge', 'wiener', etc.).
        common_support : bool
            If True, returns reconstructed densities on a shared grid.

        Examples
        --------
        >>> # 1. Setup grids and engine
        >>> x = np.linspace(-6, 6, 1000)
        >>> u = np.linspace(0, 1, 1000)
        >>> sim_engine = FDFSimulator(u_grid=u, x_grid=x)
        >>> 
        >>> # 2. Generate any custom baseline density
        >>> my_pdf = generate_base_density(x, kind='student_t', df=3)
        >>> 
        >>> # 3. Run the simulation
        >>> results = sim_engine.run_simulation(
        ...     base_pdf=my_pdf, 
        ...     n_curves=200, 
        ...     sigma=0.05, 
        ...     phis=[0.8, 0.5]
        ... )
        >>> 
        >>> # 4. Access results
        >>> simulated_densities = results['densities'] # pd.DataFrame
        """
        # 1. Forward Transform: Density -> LQD (L2 space)
        # This converts your input PDF into the Hilbert space where we add noise
        lqd_sup, lqd_base, c_val = lqdt.dens2lqd(
            dens=base_pdf, 
            dSup=self.x, 
            lqdSup=self.u, 
            t0=self.t0
        )
        
        # 2. Simulate L2 Stochastic Process
        # Generates Y_t(u) = LQD_base(u) + X_t(u) + epsilon_t(u)
        Y_l2, _, _, _, phis = simulate_l2_process(
            n=n_curves,
            base_lqd=lqd_base,
            u=self.u,
            phis=phis,
            K=dimensions,
            sigma=sigma,
            error_sigma=error_sigma,
            noise_type=noise_type,
            as_dataframe=True,
            basis=basis
        )
        self.ar_coeffs = phis
        
        # 3. Inverse Transform: LQD -> Reconstructed Densities
        # Package the simulated L2 curves with the necessary integration constants
        lqd_obj = lqdt.LQDRepresentation(
            lqd=Y_l2,
            c=c_val,
            lqd_support=self.u,
            t0=self.t0
        )
        
        reconstructed_support, df_densities = self.mlqdt.inverse_transform(
            lqd_obj, 
            common_support=common_support
        )

        self.result  = {
            "Y_l2": Y_l2,           # The simulated curves in L2 space
            "l2_mean": Y_l2.mean(axis=1),
            "densities": df_densities.fillna(0), # The final reconstructed densities
            "densities_mean": df_densities.mean(axis=1),
            "support": reconstructed_support,
            "c_val": c_val          # Useful for debugging shifts in support
        }

        return self.result
    
    def get_samples_df(self, 
                       n_samples=84
                       ) -> pd.DataFrame:
        """
        Generates an (m x T) DataFrame of samples from the simulated densities.
        m = n_samples (rows)
        T = n_curves (columns)
        
        Preserves original column names from the simulation results.

        Notes
        ---------
        * n_samples is set to 84 by default because this is the number of returns in 
        5 min intraday samples from 10:05 to 17:00 (no overnight return considered),
        just like for IBOVESPA.
        """
        if self.result is None:
            raise ValueError("No simulation results found. Run 'run_simulation' first.")
            
        df_densities = self.result["densities"]
        support = self.result["support"]
        
        # Use a dictionary comprehension for speed and column name preservation
        sampled_returns = {
            col: sim.sample_from_density(y=df_densities[col].values, x=support, n_samples=n_samples)
            for col in df_densities.columns
        }
            
        return pd.DataFrame(sampled_returns)

    def plot_simulation(self, space="density", title=None):
        if self.result is None:
            raise ValueError("No simulation results found.")

        if space == "density":
            df = self.result["densities"].copy()
            support = self.result["support"]
            y_label = "Probability Density"
            main_color = "rgba(150, 150, 150, 0.3)" 
        else:
            df = self.result["Y_l2"].copy()
            support = self.u
            y_label = "LQD Value (L2 Space)"
            main_color = "rgba(0, 100, 250, 0.3)"

        df.index = support
        
        # 1. Start with the spaghetti lines
        fig = px.line(df, template="plotly_dark")
        
        # 2. Update traces to show in legend
        # legendgroup allows you to toggle all simulated curves by clicking one
        fig.update_traces(
            line=dict(color=main_color, width=1), 
            showlegend=True,
            legendgroup="simulated" 
        )

        fig.update_layout(
            title=title or f"Simulated Functional Paths ({space.capitalize()})",
            xaxis_title="Support",
            yaxis_title=y_label,
            # hovermode="x unified",
            # Scrollable legend if the list is too long
            legend=dict(
                yanchor="top", 
                y=0.99, 
                xanchor="left", 
                x=1.02, # Moves legend outside the plot area
                traceorder="normal",
                itemsizing="constant"
            )
        )

        # Optional: Add the select menu utility you've been using
        fplt.px_select_menu(fig)
        
        return fig
    


########################################
############# GAS MODEL ################
########################################
class SkewStudentT:
    """
    Raw (non-standardized) Fernández–Steel skewed Student-t distribution.

    This class implements the asymmetric density obtained by applying
    a scale transformation to each half of a symmetric Student-t density.

    Definition
    ----------
    Let b_ν(x) be the standard Student-t density with ν degrees of freedom.
    For asymmetry parameter a > 0, define:

        h̃(x | a, ν) = c(a) * [ b_ν(a x)  if x < 0
                               b_ν(x / a) if x ≥ 0 ]

    where:
        c(a) = 2 / (a + a^{-1})

    Key properties
    --------------
    - a = 1 → symmetric Student-t
    - a > 1 → heavier right tail
    - a < 1 → heavier left tail
    - Mean ≠ 0 and variance ≠ 1 (must be standardized separately)

    This class provides:
    - density evaluation
    - raw moments required for standardization

    Example
    -------
    Evaluate the raw density and compute its moments:

    >>> import numpy as np
    >>> x = np.linspace(-4, 4, 200)
    >>> dens = SkewStudentT.density(x, a=1.5, nu=8)
    >>> mu, sigma = SkewStudentT.moments(a=1.5, nu=8)
    >>> mu, sigma
    """

    @staticmethod
    def M1_nu(nu: float) -> float:
        """
        First absolute moment over (0, ∞) of a Student-t distribution.

        This quantity appears in the analytical expression of the
        mean of the skewed-t distribution.

        Parameters
        ----------
        nu : float
            Degrees of freedom (must satisfy nu > 1)

        Returns
        -------
        float
            E[T * 1{T > 0}] where T ~ t_ν
        """
        return (np.sqrt(nu) * gamma((nu + 1) / 2)) / \
               ((nu - 1) * np.sqrt(np.pi) * gamma(nu / 2))

    @staticmethod
    def moments(a: float, nu: float):
        """
        Mean and standard deviation of the raw skewed-t distribution.

        These are used to standardize the distribution to zero mean
        and unit variance.

        Parameters
        ----------
        a : float
            Asymmetry parameter (a > 0)

        nu : float
            Degrees of freedom (nu > 2 required for variance)

        Returns
        -------
        mu : float
            Mean of the raw skewed-t

        sigma : float
            Standard deviation of the raw skewed-t
        """
        c_a = 2.0 / (a + 1.0 / a)

        # First absolute moment of symmetric t
        m1 = SkewStudentT.M1_nu(nu)

        # Variance of symmetric t
        v = nu / (nu - 2.0)

        # Analytical expressions derived from piecewise scaling
        mu = c_a * m1 * (a**2 - a**(-2))
        m2 = c_a * (v / 2.0) * (a**3 + a**(-3))

        var = m2 - mu**2
        return mu, np.sqrt(var)

    @staticmethod
    def density(x, a: float, nu: float, log: bool = False):
        """
        Evaluate the raw skewed-t density.

        Parameters
        ----------
        x : array_like
            Evaluation points

        a : float
            Asymmetry parameter

        nu : float
            Degrees of freedom

        log : bool, optional
            If True, return log-density

        Returns
        -------
        ndarray
            Density values evaluated at x
        """
        c_a = 2.0 / (a + 1.0 / a)

        # Piecewise definition implemented vectorially
        dens = np.where(
            x < 0,
            c_a * t.pdf(a * x, df=nu),
            c_a * t.pdf(x / a, df=nu)
        )

        return np.log(dens) if log else dens


class StandardizedSkewStudentT:
    """
    Standardized Fernández–Steel skewed Student-t distribution.

    This class transforms the raw skewed-t into a distribution with:
        E[X] = 0
        Var[X] = 1

    Transformation
    --------------
    If X_raw ~ h̃(x | a, ν), define:

        X_std = (X_raw - μ(a)) / σ(a)

    Then the density is:

        h(x) = σ(a) * h̃( μ(a) + σ(a) x )

    where μ(a), σ(a) come from SkewStudentT.moments.

    Notes
    -----
    - This is the object used inside the GAS likelihood.
    - Standardization ensures parameters (m, σ) have usual interpretation.

    Example
    -------
    Draw samples and evaluate density:

    >>> dist = StandardizedSkewStudentT(a=1.5, nu=8)
    >>> x = np.linspace(-3, 3, 100)
    >>> dens = dist.density(x)
    >>> sample = dist.random_numbers(1000)
    >>> np.mean(sample), np.std(sample)
    """

    def __init__(self, a, nu: float = 8.0):
        """
        Parameters
        ----------
        a : float or ndarray
            Asymmetry parameter

        nu : float
            Degrees of freedom
        """
        self.a = a
        self.nu = nu

        # Precompute standardization constants
        self.mu, self.sigma = SkewStudentT.moments(a=self.a, nu=self.nu)

    def density(self, x, log=False):
        """
        Evaluate standardized skewed-t density.

        Parameters
        ----------
        x : array_like
            Standardized input

        log : bool
            Return log-density if True

        Returns
        -------
        ndarray
        """
        # Transform back to raw scale
        y = self.mu + self.sigma * x

        logdens = (
            np.log(self.sigma) +
            SkewStudentT.density(y, a=self.a, nu=self.nu, log=True)
        )

        return logdens if log else np.exp(logdens)

    def random_numbers(self, n: int, random_state=None):
        """
        Generate i.i.d. draws from the standardized skewed-t.

        Method
        ------
        Uses the Fernández–Steel construction:

        1. Draw |T| from Student-t
        2. Assign sign with probability:
               P(X ≥ 0) = a² / (1 + a²)
        3. Apply asymmetric scaling
        4. Standardize

        Parameters
        ----------
        n : int
            Number of samples

        Returns
        -------
        ndarray
        """
        p_pos = self.a**2 / (1 + self.a**2)

        # Bernoulli draw for sign
        sign_pos = np.random.rand(n) < p_pos

        # Absolute t draws
        w = np.abs(t.rvs(df=self.nu, size=n))

        # Apply asymmetry
        x_raw = np.where(sign_pos, self.a * w, -w / self.a)

        # Standardize
        return (x_raw - self.mu) / self.sigma


class GasModel:
    """
    GAS(1,1) model with standardized skewed Student-t conditional density.

    Model structure
    ---------------
    Observation equation:
        Z_t | θ_t ~ g(· | θ_t)

    Parameterization:
        θ_t = (m_t, ℓ_t, η_t)
        σ_t = exp(ℓ_t)
        a_t = exp(η_t)

    Conditional density:
        g(z | θ) = (1/σ) * h( (z - m)/σ | a, ν )

    where h is the standardized skewed-t.

    Dynamics (γ = 0 case):
        θ_{t+1} = θ* + α ∇_t + β (θ_t - θ*)

    with:
        ∇_t = ∂ log g(Z_t | θ) / ∂θ

    Notes
    -----
    - Uses numerical differentiation for the score
    - Avoids Fisher scaling (γ = 0)
    - Suitable as a baseline computational implementation

    Example
    -------
    Full workflow: simulation + density evaluation

    >>> import numpy as np
    >>> gm = GasModel(nu=8)
    >>> sim = gm.simulate(T=200, random_state=123)
    >>> grid = np.linspace(-5, 5, 200)

    >>> dens = gm.conditional_densities(
    ...     grid=grid,
    ...     theta_path=sim["theta"]
    ... )

    >>> sim["Z"][:5]
    >>> dens.shape
    """

    def __init__(
        self,
        nu=8,
        theta_star=None,
        alpha=None,
        beta=None,
        eps: float = 1e-5,
    ):
        """
        Parameters
        ----------
        nu : float
            Degrees of freedom of skewed-t

        theta_star : ndarray, shape (3,)
            Long-run parameter level

        alpha : ndarray (3x3)
            Score sensitivity matrix

        beta : ndarray (3x3)
            Persistence matrix

        eps : float
            Step size for numerical differentiation
        """
        self.nu = nu
        self.eps = eps

        self.theta_star = (
            theta_star if theta_star is not None
            else np.array([0.0, np.log(1.0), np.log(1.0)])
        )

        self.alpha = alpha if alpha is not None else np.diag([0.03, 0.02, 0.015])
        self.beta = beta if beta is not None else np.diag([0.95, 0.90, 0.85])

        self.results = None

    def logpdf(self, z, theta):
        """
        Log conditional density log g(z | θ).

        Steps:
        1. Transform parameters to (σ, a)
        2. Standardize observation
        3. Evaluate standardized skewed-t
        4. Add Jacobian term (-log σ)

        Returns
        -------
        float or ndarray
        """
        m, ell, eta = theta

        sigma = np.exp(ell)
        a = np.exp(eta)

        x = (z - m) / sigma

        skt = StandardizedSkewStudentT(a=a, nu=self.nu)

        return -ell + skt.density(x, log=True)

    def rvs(self, n: int, theta):
        """
        Draw samples from g(z | θ).

        Implements:
            Z = m + σ * X_std

        Returns
        -------
        ndarray
        """
        m, ell, eta = theta

        sigma = np.exp(ell)
        a = np.exp(eta)

        skt = StandardizedSkewStudentT(a=a, nu=self.nu)

        return m + sigma * skt.random_numbers(n)

    def score(self, z, theta):
        """
        Numerical score vector ∇_t.

        Uses central finite differences:
            ∂f/∂θ ≈ [f(θ+h) - f(θ-h)] / (2h)

        Adaptive step size improves stability.

        Returns
        -------
        ndarray, shape (3,)
        """
        theta = np.asarray(theta)
        grad = np.zeros_like(theta)

        for j in range(len(theta)):
            step = self.eps * (1 + abs(theta[j]))

            theta_plus = theta.copy()
            theta_minus = theta.copy()

            theta_plus[j] += step
            theta_minus[j] -= step

            lp = self.logpdf(z, theta_plus)
            lm = self.logpdf(z, theta_minus)

            grad[j] = (lp - lm) / (2 * step)

        return grad

    def simulate(
            self, 
            T: int, 
            theta_init = None, 
            burn_in : int = 100, 
            random_state = None
            ):
        """
        Simulate a trajectory from the GAS model.

        Algorithm
        ---------
        For t = 1,...,T:
            1. Draw Z_t | θ_t
            2. Compute score ∇_t
            3. Update θ_{t+1}

        Returns
        -------
        dict containing:
            Z           : observations
            theta       : parameter path
            mean        : m_t
            sd          : σ_t
            asymmetry   : a_t
            score       : score sequence
        """
        if random_state is not None:
            np.random.seed(random_state)
        if burn_in:
            T += burn_in

        theta_init = (
            theta_init.copy()
            if theta_init is not None
            else self.theta_star.copy()
        )

        k = len(self.theta_star)

        Z = np.zeros(T)
        theta_path = np.zeros((T + 1, k))
        scores = np.zeros((T, k))

        theta_path[0] = theta_init

        for t in range(T):
            theta_t = theta_path[t]

            # Step 1: draw observation
            Z[t] = self.rvs(1, theta_t)[0]

            # Step 2: compute score
            grad = self.score(Z[t], theta_t)
            scores[t] = grad

            # Step 3: GAS update
            theta_path[t + 1] = (
                self.theta_star
                + self.alpha @ grad
                + self.beta @ (theta_t - self.theta_star)
            )

        results = {
                    "Z": Z[burn_in:],
                    "theta": theta_path[burn_in:-1],
                    "mean": theta_path[burn_in:-1, 0],
                    "sd": np.exp(theta_path[burn_in:-1, 1]),
                    "asymmetry": np.exp(theta_path[burn_in:-1, 2]),
                    "score": scores[burn_in:],
                }
        
        self.results = results
        return results

    def conditional_densities(self, grid, theta_path, log=False):
        """
        Evaluate f_t(z) over a grid for all t.

        Fully vectorized over:
            - time (rows)
            - grid (columns)

        Useful for:
        - visualization
        - density forecast evaluation

        Returns
        -------
        DataFrame
            Rows = grid
            Columns = time index
        """
        grid = np.asarray(grid)
        theta_path = np.asarray(theta_path)

        m = theta_path[:, 0][:, None]
        ell = theta_path[:, 1][:, None]
        eta = theta_path[:, 2][:, None]

        sigma = np.exp(ell)
        a = np.exp(eta)

        x = (grid[None, :] - m) / sigma

        skt = StandardizedSkewStudentT(a=a, nu=self.nu)
        logdens = -ell + skt.density(x, log=True)

        return pd.DataFrame(
            logdens if log else np.exp(logdens)
        ).T.set_index(grid)
        
    def bootstrap(self, sample: pd.Series, grid: np.ndarray, n_boot: int = 1000):
            """
            Executes the bootstrap procedure for time T.
            
            Parameters
            ----------
            sample : np.ndarray
                The cross-sectional data points (e.g., 288 points) observed 
                at the final time step T-1.
            n_boot : int
                Number of bootstrap iterations (B).

            Returns
            -------
            list of dicts
                Each dict contains 'theta_t_b' (the bootstrap density parameters)
                and 'Z_it_b' (the simulated cross-section for that density).
            """
            date = sample.name
            sample = sample.values

            if self.results is None:
                raise ValueError("You must run .simulate() before bootstrapping.")

            # Step 0: Extract theta_{T-1} from the last state of simulation
            # theta_path is (T, k), so [-1] is the final parameter state
            theta_prev = self.results["theta"][-1]
            n_t = len(sample)
            
            boot_results = []

            for b in range(n_boot):
                # Step (a): Randomly sample one z from the cross-section Z_{t-1}
                # and generate the bootstrap density parameters theta_t^(b)
                z_sample = np.random.choice(sample)
                
                # Compute score using the sampled observation and the anchor theta
                grad = self.score(z_sample, theta_prev)
                
                # Update to obtain f_t^(b) parameters
                theta_t_b = (
                    self.theta_star
                    + self.alpha @ grad
                    + self.beta @ (theta_prev - self.theta_star)
                )

                # Step (b): Sample n_t points from the resulting density f_t^(b)
                Z_it_b = self.rvs(n_t, theta_t_b)
                Z_it_b = pd.DataFrame(Z_it_b)
                Z_it_b.columns = [date]

                density_df = self.conditional_densities(
                                grid=grid, 
                                theta_path=theta_t_b.reshape(1, -1)
                            )

                boot_results.append({
                    "theta_t_b": theta_t_b,
                    "density": density_df,
                    "Z_it_b": Z_it_b
                })

            return boot_results
    