import numpy as np
import pandas as pd
from typing import Tuple

def rv(
    df_returns: pd.DataFrame,
    log_rv: bool = False,
    name:   str  = None
    ) -> Tuple[pd.Series, pd.DataFrame]:
    """
    Compute realized variance (RV) from sampled intraday returns.

    The function assumes that each column in `df_returns` corresponds to a time
    period (e.g., a day), and each row represents an intraday return (e.g.,
    5-minute returns) sampled from a forecasted density.

    Realized variance is computed as:

        RV_t = sum_{i=1}^m r_{i,t}^2

    where:
        - m is the number of intraday observations
        - r_{i,t} is the i-th intraday return at time t

    Parameters
    ----------
    df_returns : pd.DataFrame
        DataFrame of sampled returns with shape (m, T), where:
        - m : number of intraday samples (rows)
        - T : number of time periods (columns)
        Columns represent time indices (e.g., dates).

    log_rv : bool, default=False
        If True, returns the logarithm of realized variance.

    Returns
    -------
    Tuple[pd.Series, pd.DataFrame]
        realized_variance : pd.Series
            Realized variance for each time period (length T).
        df_squared_returns : pd.DataFrame
            Squared intraday returns with same shape as input.
    """

    if df_returns is None or df_returns.empty:
        raise ValueError("Input DataFrame 'df_returns' must not be empty.")

    # Ensure numeric data
    if not np.issubdtype(df_returns.values.dtype, np.number):
        raise TypeError("Input DataFrame must contain numeric values.")

    # Squared returns
    df_squared_returns: pd.DataFrame = df_returns.pow(2)

    # Realized variance (sum over intraday dimension)
    realized_variance: pd.Series = df_squared_returns.sum(axis=0)

    if log_rv:
        realized_variance = np.log(realized_variance)

    # Naming for clarity in downstream pipelines
    if name is not None:
        realized_variance.name = name
    else:
        realized_variance.name = "realized_volatility"
    realized_variance.index.name = "date"
    
    return realized_variance, df_squared_returns