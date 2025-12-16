import matplotlib.pyplot as plt
import pandas as pd
from typing import Dict
import numpy as np

from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import adfuller, kpss#, phillips_perron

def adf_test(x):
    res = adfuller(x, autolag='AIC')
    # print('ADF Statistic: %f' % res[0])
    # print('p-value: %f' % res[1])
    # print('Critical Values:')
    # for k,v in res[4].items():
    #     print(f'   {k}: {v}')
    if not res[1] < 0.05:
        print('\t ADF => reject H0 (unit root) ? ', res[1] < 0.05)
        plt.figure()
        plt.plot(x)
        plt.show()
    print()

def kpss_test(x, regression='c'):
    statistic, p_value, lags, crit = kpss(x, regression=regression)
    # print('KPSS Statistic: %f' % statistic)
    # print('p-value: %f' % p_value)
    # print('Critical Values:')
    # for k,v in crit.items():
    #     print(f'   {k}: {v}')
    if p_value < 0.05:
        print('\t KPSS => reject H0 (stationary) ? ', p_value < 0.05)
        plt.figure()
        plt.plot(x)
        plt.show()
    print()

    
def select_order_ic(
        data: pd.DataFrame, 
        maxlags: int = 10
        ) -> Dict[str, int]:
    """
    Select lag order using AIC and BIC from statsmodels VAR.select_order.
    Returns dict with keys 'aic', 'bic', 'hqic' (if available).
    """
    model = VAR(data)
    sel = model.select_order(maxlags)
    # statsmodels returns object with attributes aic, bic, hqic that are integers (lags)
    return {'aic': int(sel.aic), 'bic': int(sel.bic), 'hqic': int(sel.hqic)}

def fit_var(data: pd.DataFrame, nlags: int) -> VAR:
    """
    Fit a VAR model and return the fitted results object.
    """
    model = VAR(data)
    res = model.fit(nlags)
    return res

def forecast_var(res, steps: int = 1) -> np.ndarray:
    """
    Forecast using a fitted statsmodels VARResults object.
    Returns a numpy array (steps x k).
    """
    return res.forecast(res.endog[-res.k_ar:], steps=steps)









############################################################################################
##################################### ACCURACY METRICS #####################################
############################################################################################
 
def KLdiv_matrix(object, eps=1e-4, overlap=True):
    """
    Kullback–Leibler divergence between columns of a matrix.
    Adapted from R package flexmix and function KLdiv.

    Parameters
    ----------
    object : np.ndarray, shape (n, p)
        Each column is a discrete density / pmf.
    eps : float
        Small positive constant to avoid log(0).
    overlap : bool
        If False, requires positive overlap between densities.

    Returns
    -------
    z : np.ndarray, shape (p, p)
        Pairwise KL divergence matrix.
    """
    if not np.issubdtype(object.dtype, np.number):
        raise ValueError("object must be a numeric matrix")

    object = object.copy()
    p = object.shape[1]

    # initialize result
    z = np.full((p, p), np.nan)

    # replace small values
    w = object < eps
    if np.any(w):
        object[w] = eps

    # normalize columns (like sweep(..., 2, colSums, "/"))
    object = object / object.sum(axis=0, keepdims=True)

    for k in range(p - 1):
        for l in range(1, p):
            ok = (object[:, k] > eps) & (object[:, l] > eps)
            if (not overlap) or np.any(ok):
                z[k, l] = np.sum(
                    object[:, k] *
                    (np.log(object[:, k]) - np.log(object[:, l]))
                )
                z[l, k] = np.sum(
                    object[:, l] *
                    (np.log(object[:, l]) - np.log(object[:, k]))
                )

    np.fill_diagonal(z, 0.0)
    return z

def KLdiv(
        p : np.array, 
        q : np.array, 
        eps=1e-4, 
        overlap=True) -> np.float64:
    """
    Kullback–Leibler divergence KL(p || q) for two discrete densities.
    Adapted from R package flexmix and function KLdiv.

    Parameters
    ----------
    p, q : array-like, shape (n,)
        Discrete densities defined on the same grid.
    eps : float
        Small positive constant to avoid log(0).
    overlap : bool
        If False, requires positive overlap between p and q
        (same logic as the original R code).

    Returns
    -------
    float
        KL(p || q)

    Example
    -------
    x = np.linspace(-3, 3, 200)
    from scipy.stats import norm, uniform, t
    y = pd.DataFrame(np.array([uniform.pdf(x), norm.pdf(x), t.pdf(x, df=10)]).T, columns=["Uniform", "Normal", "t"])
    y.plot()
    KLdiv(y.Uniform, y.Normal)
    KLdiv(y.t, y.Normal) # more similar than Uniform and Normal
    """
    p = np.asarray(p, dtype=float).copy()
    q = np.asarray(q, dtype=float).copy()

    if p.ndim != 1 or q.ndim != 1:
        raise ValueError("p and q must be one-dimensional arrays")
    if p.shape[0] != q.shape[0]:
        raise ValueError("p and q must have the same length")

    # replace small values
    p[p < eps] = eps
    q[q < eps] = eps

    # renormalize (same as sweep + colSums)
    p /= p.sum()
    q /= q.sum()

    # overlap condition (faithful to the R code)
    ok = (p > eps) & (q > eps)
    if overlap and not np.any(ok):
        return np.nan
    
    kld = np.sum(p * (np.log(p) - np.log(q)))

    return kld

def JSdiv(
        p: np.array,
        q: np.array,
        eps: float = 1e-4,
        overlap: bool = True
    ) -> np.float64:
    """
    Jensen–Shannon divergence between two discrete probability densities.
    Adapted from Kokoszka 2019.

    This function computes the (symmetric) Jensen–Shannon divergence (JSD)
    between two discrete densities `p` and `q` defined on a common grid.
    The JSD is defined as

        JSD(p, q) = 0.5 * KL(p || m) + 0.5 * KL(q || m),

    where

        m = 0.5 * (p + q)

    is the arithmetic mean mixture of the two densities, and `KL(·||·)`
    denotes the Kullback–Leibler divergence.

    The implementation follows standard practice in density forecast
    evaluation and information theory, using thresholding and
    renormalization to ensure numerical stability.

    Parameters
    ----------
    p, q : array-like of shape (n,)
        Discrete probability densities evaluated on the same support/grid.
        The inputs need not be perfectly normalized; normalization is
        enforced internally.

    eps : float, optional (default=1e-4)
        Small positive constant used to threshold the densities in order
        to avoid undefined logarithms. Values smaller than `eps` are
        replaced by `eps` prior to normalization.

    overlap : bool, optional (default=True)
        If True, the overlap condition is enforced inside the KL divergence
        computations, mimicking the behavior of the original R
        implementation from the `flexmix` package. If the overlap condition
        is violated, the function returns `np.nan`.

    Returns
    -------
    jsd : np.float64
        Jensen–Shannon divergence between `p` and `q`. The value is
        non-negative and finite. Using natural logarithms, the divergence
        is bounded above by `log(2)`.

    Notes
    -----
    - The Jensen–Shannon divergence is symmetric, i.e.,
      JSD(p, q) = JSD(q, p).
    - Unlike the Kullback–Leibler divergence, JSD is always finite when
      computed with the arithmetic mean mixture.
    - The square root of JSD defines a metric on the space of probability
      distributions.
    - This implementation relies on the `KLdiv` function and inherits its
      numerical and overlap-handling conventions.

    See Also
    --------
    KLdiv : Kullback–Leibler divergence for discrete densities.

    Examples
    --------
    x = np.linspace(-3, 3, 200)
    from scipy.stats import norm, uniform, t
    y = pd.DataFrame(np.array([uniform.pdf(x), norm.pdf(x), t.pdf(x, df=10)]).T, columns=["Uniform", "Normal", "t"])
    y.plot()
    JSdiv(y.Uniform, y.Normal) 
    JSdiv(y.t, y.Normal) # more similar than Uniform and Normal
    """
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)

    # Threshold to avoid log(0)
    p[p < eps] = eps
    q[q < eps] = eps

    # Renormalize to ensure valid probability densities
    p /= p.sum()
    q /= q.sum()

    # Arithmetic mean mixture
    m = 0.5 * (p + q)

    jsd = 0.5 * KLdiv(p, m, eps=eps, overlap=overlap) \
        + 0.5 * KLdiv(q, m, eps=eps, overlap=overlap)

    return jsd

def L_norm(
        p: np.array,
        q: np.array,
        norm_type: str = "L1"
    ) -> np.float64:
    """
    Discrete Lp-type norm between two functions evaluated on a common grid.

    Parameters
    ----------
    p, q : array-like of shape (n,)
        Discrete evaluations of two functions (e.g., densities) on the
        same support/grid.

    norm_type : {"L1", "L2", "L_inf"}, optional (default="L1")
        Type of norm to compute:
        - "L1"    : sum_i |p_i - q_i|
        - "L2"    : sqrt( sum_i (p_i - q_i)^2 )
        - "LINF" : max_i |p_i - q_i|

    Returns
    -------
    norm : np.float64
        Value of the selected discrete norm.

    Notes
    -----
    - These are discrete norms; no integration weights are applied.
    - The norms treat densities as generic functions and do not exploit
      their probabilistic structure.
    """
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)

    if p.shape != q.shape:
        raise ValueError("p and q must have the same shape")

    if norm_type == "L1":
        norm = np.sum(np.abs(p - q))

    elif norm_type == "L2":
        norm = np.sqrt(np.sum((p - q) ** 2))

    elif norm_type == "LINF":
        norm = np.max(np.abs(p - q))

    else:
        raise ValueError("norm_type must be one of {'L1', 'L2', 'LINF'}")
    
    return norm
    

def get_metrics(
        p: np.array,
        q: np.array
    ) -> Dict:
    """
    Compute a collection of functional and divergence-based error measures
    between two discrete densities.

    This function evaluates multiple discrepancy measures between two
    density-valued functions `p` and `q` defined on a common grid. The
    returned metrics include both information-theoretic divergences and
    classical functional norms, allowing for a comprehensive assessment
    of density forecast accuracy.

    The following measures are computed:
    - Kullback–Leibler divergence (KLD)
    - Jensen–Shannon divergence (JSD)
    - L1 norm (integrated absolute error, discrete)
    - L2 norm (root integrated squared error, discrete)
    - L-infinity norm (supremum norm, discrete)

    Parameters
    ----------
    p, q : array-like of shape (n,)
        Discrete probability densities evaluated on the same support/grid.
        The inputs need not be perfectly normalized; normalization and
        numerical stabilization are handled internally by the divergence
        functions.

    Returns
    -------
    metrics_dict : dict
        Dictionary containing the following key–value pairs:

        - "KLD"       : float
            Kullback–Leibler divergence KL(p || q).
        - "JSD"       : float
            Jensen–Shannon divergence between p and q.
        - "L1_norm"   : float
            Discrete L1 norm ||p − q||₁.
        - "L2_norm"   : float
            Discrete L2 norm ||p − q||₂.
        - "LINF_norm" : float
            Discrete L-infinity norm ||p − q||∞.

    Notes
    -----
    - KLD and JSD explicitly exploit the probabilistic structure of
      densities and quantify distributional discrepancy and information
      loss.
    - L1, L2, and L-infinity norms treat densities as generic functions
      and measure pointwise approximation error on the grid.
    - All norms are discrete and depend on the grid resolution; no
      integration weights are applied.
    - This function is intended for forecast evaluation in functional
      time series of densities, where complementary error measures
      provide a more complete comparison than any single metric.

    See Also
    --------
    KLdiv : Kullback–Leibler divergence for discrete densities.
    JSdiv : Jensen–Shannon divergence.
    L_norm : Dispatcher for discrete Lp-type norms.

    Examples
    --------
    x = np.linspace(-3, 3, 200)
    from scipy.stats import norm, uniform, t
    y = pd.DataFrame(np.array([uniform.pdf(x), norm.pdf(x), t.pdf(x, df=10)]).T, columns=["Uniform", "Normal", "t"])
    metrics = get_metrics(y.Uniform, y.Normal)
    metrics["JSD"]
    metrics["L2_norm"]
    """
    KLD  = KLdiv(p, q)
    JSD  = JSdiv(p, q)
    L1   = L_norm(p, q, norm_type="L1")
    L2   = L_norm(p, q, norm_type="L2")
    LINF = L_norm(p, q, norm_type="LINF")

    metrics_dict = {
        "KLD": KLD,
        "JSD": JSD,
        "L1_norm":   L1,
        "L2_norm":   L2,
        "LINF_norm": LINF
    }

    return metrics_dict

def overall_KLD(
        X : pd.DataFrame,
        Y : pd.DataFrame,
        n_test : int
    ):
    """
    Compute the overall symmetric Kullback–Leibler divergence between
    two sequences of discretized density functions.

    For each evaluation index, this function computes the Kullback–Leibler
    divergence between the observed density and its forecast in both
    directions, KL(f || f̂) and KL(f̂ || f). The resulting divergences are
    then averaged over the evaluation sample and summed, yielding the
    time-averaged Jeffreys divergence as an overall measure of density
    forecast accuracy.

    Parameters
    ----------
    X : pandas.DataFrame
        Observed densities evaluated on a common grid. Each column
        corresponds to a time index (or forecast origin), and each
        column is assumed to represent a discrete probability density
        function.

    Y : pandas.DataFrame
        Forecast densities evaluated on the same grids as `X`.
        Must have the same shape and column ordering as `X`.

    n_test : int
        Number of evaluation points (columns) over which the divergence
        is computed.

    Returns
    -------
    summary : float
        Overall density forecast error measured as the average symmetric
        Kullback–Leibler divergence (Jeffreys divergence) across the
        evaluation period.

    Notes
    -----
    * Divergences are computed independently for each column (time index).
    * The function `KLdiv` is assumed to compute KL(p || q) for two
      discrete densities p and q.
    * The final summary is given by
        (1 / n_test) * sum_t [ KL(f_t || f̂_t) + KL(f̂_t || f_t) ],
      where f_t and f̂_t denote the observed and forecast densities at
      time t.
    * The resulting measure is symmetric but is not itself a proper
      scoring rule.
    """
    
    measures = np.full((n_test, 2), np.nan)
    for p in range(n_test):
        p_true = X.iloc[:, p]
        p_fc = Y.iloc[:, p]
        measures[p, 0] = KLdiv(p_true, p_fc)
        measures[p, 1] = KLdiv(p_fc, p_true)

    summary = np.round(np.nanmean(measures, axis=0).sum(), 4)

    return summary

def overall_JSD(
        X: pd.DataFrame,
        Y: pd.DataFrame,
        n_test: int
    ) -> float:
    """
    Compute the overall Jensen–Shannon divergence between two sequences
    of discretized density functions using the arithmetic mean mixture.

    For each evaluation index, the Jensen–Shannon divergence between the
    observed and forecast densities is computed. The resulting values
    are then summed over the evaluation sample to obtain an overall
    measure of density forecast accuracy.

    Parameters
    ----------
    X : pandas.DataFrame
        Observed densities evaluated on a common grid. Each column
        corresponds to a time index.

    Y : pandas.DataFrame
        Forecast densities evaluated on the same grid as `X`.
        Must have the same shape and column ordering as `X`.

    n_test : int
        Number of evaluation points (columns).

    Returns
    -------
    summary : float
        Overall Jensen–Shannon divergence summed over the evaluation
        period.

    Notes
    -----
    * The Jensen–Shannon divergence is computed as
          JSD(f_t, f̂_t) =
          0.5 * KL(f_t || m_t) + 0.5 * KL(f̂_t || m_t),
      where m_t = 0.5 * (f_t + f̂_t).
    * Divergences are computed independently for each time index and
      aggregated by summation, matching the original R implementation.
    * The resulting measure is symmetric, finite, and bounded above
      by log(2) when natural logarithms are used.

    See Also
    --------
    KLdiv : Kullback–Leibler divergence for discrete densities.
    JSdiv : Jensen–Shannon divergence for discrete densities.
    """

    jsd = np.full(n_test, np.nan)

    for t in range(n_test):
        p_true = X.iloc[:, t].to_numpy()
        p_fc   = Y.iloc[:, t].to_numpy()

        jsd[t] = JSdiv(p_true, p_fc)

    summary = np.round(np.nansum(jsd), 4)

    return summary

def overall_Lnorm(
        X: pd.DataFrame,
        Y: pd.DataFrame,
        n_test: int,
        norm: str = "L1"
    ) -> float:
    """
    Compute the overall L^p-type norm between observed and forecast
    density functions over an evaluation sample.

    Parameters
    ----------
    X : pandas.DataFrame
        Observed densities evaluated on a common grid.
        Each column corresponds to a time index.

    Y : pandas.DataFrame
        Forecast densities evaluated on the same grid as `X`.

    n_test : int
        Number of evaluation points (columns).

    norm : {"L1", "L2", "LINF"}, optional
        Type of norm to compute:
        - "L1"   : sum of absolute deviations
        - "L2"   : square root of sum of squared deviations
        - "LINF" : maximum absolute deviation

    Returns
    -------
    summary : float
        Mean norm value over the evaluation period, rounded to 4 decimals.

    Notes
    -----
    - Norms are computed pointwise for each time index and then averaged
      over the evaluation sample.
    - This function reproduces the aggregation used in the original R
      implementation exactly.
    """

    values = np.full(n_test, np.nan)

    for t in range(n_test):
        diff = X.iloc[:, t].to_numpy() - Y.iloc[:, t].to_numpy()

        if norm.upper() == "L1":
            values[t] = np.sum(np.abs(diff))

        elif norm.upper() == "L2":
            values[t] = np.sqrt(np.sum(diff ** 2))

        elif norm.upper() == "LINF":
            values[t] = np.max(np.abs(diff))

        else:
            raise ValueError("norm must be one of {'L1', 'L2', 'LINF'}")

    return np.round(np.nanmean(values), 4)