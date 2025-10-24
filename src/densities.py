import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

filepath = r"C:\Users\user\Projetos\densities4risk\data\processed\database.parquet"
df = pd.read_parquet(filepath)

TICKER = 'BVSP'

df_bvsp = df[df["ticker"]==TICKER].loc[:,["datetime", "open", "high", "low", "close"]].dropna().drop_duplicates()
df_bvsp.set_index("datetime", inplace=True)

start_time = '10:00'
end_time = '17:00'
freq = '5T'  # 5 minutes

# Create a complete time index for each day
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

# Interpolar valores vazios linearmente
df_complete.interpolate(method="linear", inplace=True)

# Retornos
df_complete["R_t"] = df_complete["close"].pct_change()
df_complete["r_t"] = np.log(df_complete["close"]) - np.log(df_complete["close"].shift(1))

df_complete = df_complete.dropna(subset="r_t")

# Group by day
grouped = df_complete.groupby(df_complete['datetime'].dt.date)

# Compute KDE per day and store results
# kde_results = {}
densities = {}

for day, group in grouped:
    data = group['r_t'].dropna()
    if len(data) > 1:  # need at least 2 points
        kde = gaussian_kde(data)
        # Evaluate on a grid for visualization
        x_grid = np.linspace(data.min(), data.max(), 200)
        y_kde = kde(x_grid)
        # kde_results[day] = (x_grid, y_kde)
        densities[day]=y_kde

df_densities = pd.DataFrame(densities)

df_densities.to_excel("data/processed/kde.xlsx", index=False)