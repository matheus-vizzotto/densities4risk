import pandas as pd
import numpy as np


# filepath = r"C:\Users\user\Projetos\densities4risk\data\processed\database.parquet"
# df = pd.read_parquet(filepath)


# TICKER = 'BVSP'
# df_bvsp = df[df["ticker"]==TICKER].loc[:,["datetime", "open", "high", "low", "close"]].dropna().drop_duplicates()
# df_bvsp.set_index("datetime", inplace=True)


# # Create a complete time index for each day
# start_time = '10:00'
# end_time = '17:00'
# freq = '5T'  # 5 minutes
# def complete_day(group):
#     day = group.index[0].date()
#     full_index = pd.date_range(
#         f"{day} {start_time}",
#         f"{day} {end_time}",
#         freq=freq
#     )
#     return group.reindex(full_index)

# # Apply to each day
# df_complete = df_bvsp.groupby(df_bvsp.index.date, group_keys=False).apply(complete_day)
# # Reset index and rename column
# df_complete = df_complete.reset_index().rename(columns={'index': 'datetime'})
# df_complete.interpolate(method="linear", inplace=True)


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