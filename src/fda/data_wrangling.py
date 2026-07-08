import numpy as np
import pandas as pd

def complete_day(group, start_time='10:00', end_time='17:00', freq='5T'):
    """Create a complete time index for each day
    """
    day = group.index[0].date()
    full_index = pd.date_range(
        f"{day} {start_time}",
        f"{day} {end_time}",
        freq=freq
    )
    return group.reindex(full_index)


def ts_to_df(
        df : pd.DataFrame, 
        y : str, 
        datetime_col : str
        ) -> pd.DataFrame:
    """Transforms a time series dataframe to a wide format where columns are dates and
    rows are times. 

    Args:
        df (pd.DataFrame): _description_
        y (str): _description_
        datetime_col (str): _description_

    Returns:
        pd.DataFrame: _description_
    """
    df["date"] = df[datetime_col].dt.date
    df["time"] = df[datetime_col].dt.time
    df_wide = df.pivot(index="time", columns="date", values= y)
    df_wide.fillna(df_wide.mean(axis=0).mean(), inplace=True)

    return df_wide

def add_return_cols(
        df: pd.DataFrame, 
        y:  str,
        dropna: bool=True
        ) -> pd.DataFrame:
    """
    Calculate simple and logarithmic returns for a financial time series.
    
    This function computes both simple percentage returns (R_t) and continuous 
    logarithmic returns (r_t) for a specified price column in a DataFrame. 
    Both return types are added as new columns to a copy of the original DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing the financial time series data.
        Must contain the column specified by `y`.
    y : str
        Name of the column containing price/asset values for which returns 
        should be calculated. This column should typically contain closing 
        prices or similar financial metrics.
    dropna : bool, default=True
        If True, removes rows with NaN values (typically the first row after 
        calculating returns since it has no previous value for differencing).
        If False, retains all rows with NaN in the return columns.
    
    Returns
    -------
    pd.DataFrame
        A new DataFrame containing all original columns plus two new columns:
        - 'R_t': Simple percentage returns calculated as pct_change()
        - 'r_t': Logarithmic returns calculated as log(P_t) - log(P_{t-1})
    
    Notes
    -----
    - Simple returns (R_t): (P_t - P_{t-1}) / P_{t-1}
    - Logarithmic returns (r_t): ln(P_t) - ln(P_{t-1}) ≈ ln(1 + R_t)
    """
    
    R_t = df.loc[:,y].pct_change()
    r_t = np.log(df.loc[:, y]) - np.log(df.loc[:,y].shift(1))
    
    df2 = df.copy(deep=True)
    df2["R_t"] = R_t
    df2["r_t"] = r_t

    if dropna:
        df2.dropna(inplace=True)

    return df2


def build_lags(
        df : pd.DataFrame,
        col : str,
        n_lags : int,
        dropna : bool = True
    ) -> pd.DataFrame:

    lags = []
    for lag in range(n_lags+1):
        lags.append(df.loc[:, col].shift(lag))

    df = pd.concat(lags, axis=1)
    df_columns = [f"t-{i}" for i in range(n_lags+1)]
    df.columns = df_columns

    if dropna:
        df.dropna(inplace=True)

    return df

def mark_repeated_columns(df, threshold=10, dropna=True):
    """
    Identify columns that contain any value repeated more than `threshold` times.
    Function to detect weird data, like returns on Ash Wednesday (started on 13pm).
    Already detected: '2025-03-05', '2025-08-01'.

    Parameters
    ----------
    df : pd.DataFrame
        Columns represent days.
    threshold : int
        Maximum allowed repetitions of a single value.
    dropna : bool
        Whether to ignore NaNs when counting repetitions.

    Returns
    -------
    bad_cols : list
        Column names with excessive repetition.
    mask : pd.Series
        Boolean mask indexed by columns.
    """

    def has_too_many_repeats(col):
        counts = col.value_counts(dropna=dropna)
        return counts.max() > threshold if not counts.empty else False

    mask = df.apply(has_too_many_repeats, axis=0)

    bad_cols = mask[mask].index.tolist()

    return bad_cols, mask