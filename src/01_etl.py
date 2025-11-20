import pandas as pd
import numpy as np

database = pd.read_parquet('data/interim/database.parquet')

# Série específica
TICKER = 'BVSP'
df_bvsp = database[database["ticker"]==TICKER].loc[:,["datetime", "open", "high", "low", "close"]].dropna().drop_duplicates()
df_bvsp.set_index("datetime", inplace=True)
start_time = '10:00'
end_time = '17:00'
freq = '5T'  # 5 minutes

def complete_day(group):
    """Create a complete time index for each day
    """
    day = group.index[0].date()
    full_index = pd.date_range(
        f"{day} {start_time}",
        f"{day} {end_time}",
        freq=freq
    )
    return group.reindex(full_index)

df_complete = df_bvsp.groupby(df_bvsp.index.date, group_keys=False).apply(complete_day)
df_complete = df_complete.reset_index().rename(columns={'index': 'datetime'})
df_complete.interpolate(method="linear", inplace=True)
df_complete["R_t"] = df_complete["close"].pct_change()
df_complete["r_t"] = np.log(df_complete["close"]) - np.log(df_complete["close"].shift(1))
df_complete.dropna(inplace=True)

df_complete.to_excel("data/processed/BVSP_returns.xlsx", index=False)

# Saves file in functional data format
df = df_complete.loc[:,["datetime", "r_t"]]
df["date"] = df["datetime"].dt.date
df["time"] = df["datetime"].dt.time
returns = df.pivot(index="time", columns="date", values="r_t")
returns.fillna(returns.mean(axis=0).mean(), inplace=True)
returns.to_excel("data/processed/BVSP_returns_wide.xlsx")