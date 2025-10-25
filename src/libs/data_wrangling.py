import pandas as pd


filepath = r"C:\Users\user\Projetos\densities4risk\data\processed\database.parquet"
df = pd.read_parquet(filepath)


TICKER = 'BVSP'
df_bvsp = df[df["ticker"]==TICKER].loc[:,["datetime", "open", "high", "low", "close"]].dropna().drop_duplicates()
df_bvsp.set_index("datetime", inplace=True)


# Create a complete time index for each day
start_time = '10:00'
end_time = '17:00'
freq = '5T'  # 5 minutes
def complete_day(group):
    day = group.index[0].date()
    full_index = pd.date_range(
        f"{day} {start_time}",
        f"{day} {end_time}",
        freq=freq
    )
    return group.reindex(full_index)

# Apply to each day
df_complete = df_bvsp.groupby(df_bvsp.index.date, group_keys=False).apply(complete_day)
# Reset index and rename column
df_complete = df_complete.reset_index().rename(columns={'index': 'datetime'})
df_complete.interpolate(method="linear", inplace=True)