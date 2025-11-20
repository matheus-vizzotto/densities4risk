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

def calc_returns(
        df : pd.DataFrame, 
        y : str) -> pd.DataFrame:
    """Adds percantual and log returns columns to the original dataframe.

    Args:
        df (pd.DataFrame): _description_
        y (str): _description_

    Returns:
        pd.DataFrame: _description_
    """

    df["R_t"] = df[y].pct_change()
    df["r_t"] = np.log(df["close"]) - np.log(df["close"].shift(1))
    df.dropna(inplace=True)

    return df