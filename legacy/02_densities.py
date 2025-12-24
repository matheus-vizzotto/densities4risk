import numpy as np
import pandas as pd
import json

from scipy.stats import gaussian_kde

import matplotlib.pyplot as plt
import plotly.graph_objects as go

df = pd.read_excel("data/processed/BVSP_returns.xlsx")

# Group by day
grouped = df.groupby(df['datetime'].dt.date)

# Compute KDE per day and store results
kdes = {}
day_index = 1
r_minus_total = []
r_plus_total = []
for day, group in grouped:
    # y variable
    data = group['r_t']
    # kde object
    kde = gaussian_kde(data)
    # support
    r_minus = data.min()
    r_plus = data.max()
    r_minus_total.append(r_minus)
    r_plus_total.append(r_plus)
    # grid
    x_grid = np.linspace(data.min(), data.max(), 200)
    # kde estimation
    y_kde = kde(x_grid)
    # store values
    kdes[day_index] = {
        "x_grid": x_grid.tolist(),
        "y_kde":  y_kde.tolist()
        }
    
    day_index += 1

file_path = "data/processed/densities.json"
with open(file_path, 'w') as json_file:
    json.dump(kdes, json_file, indent=4)




# 1) Global support
df_returns = pd.read_excel("data/processed/BVSP_returns_wide.xlsx", index_col="time")

global_min = df_returns.min().min()
global_max = df_returns.max().max()
m = 3000
u = np.linspace(global_min, global_max, m)

# 2) Prepare density matrix (m × T)
df_densities = pd.DataFrame(index=u, columns=df_returns.columns)

# 3) KDE for each day evaluated on a common support
for t in df_returns.columns:
    kde = gaussian_kde(df_returns[t])
    df_densities[t] = kde(u)

df_densities.to_excel("data/processed/BVSP_returns_densities.xlsx")


############################################# VISUALIZAÇÃO ######################################
# Distribuição de retornos
plt.figure(figsize=(15,5))
plt.hist(df["r_t"], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
plt.title("Distribution of log-returns")
plt.savefig("img/distribution_log_returns.png")

# Endpoints dos suportes de densidades diárias
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.hist(r_minus_total, bins=30, alpha=0.7, color='skyblue', edgecolor='black', label="$r^-$")
ax1.legend()
ax2.hist(r_plus_total, bins=30, alpha=0.7, color='lightcoral', edgecolor='black', label="$r^+$")
ax2.legend()
plt.suptitle('Densities supports endpoints')
plt.tight_layout()
plt.savefig('img/densidies_supports.png')

# Densidades estimadas sobrepostas
plt.figure(figsize=(15,5))
for day_index in kdes.keys():
    plt.plot(kdes[day_index]["x_grid"], kdes[day_index]["y_kde"], linewidth=1, alpha=.2, c = '#8338ec')
plt.title("Kernel Density Estimates")
plt.savefig('img/kdes.png')

# Densidades estimadas 3d
days_sorted = sorted(kdes.keys())
all_x = np.concatenate([kdes[d]["x_grid"] for d in days_sorted])
x_common = np.linspace(all_x.min(), all_x.max(), 200)
Z = []
for d in days_sorted:
    kde = kdes[d]["x_grid"]
    density = kdes[d]["y_kde"]
    # Interpolate densities to common x grid
    z_interp = np.interp(x_common, kde, density)
    Z.append(z_interp)
Z = np.array(Z)
Y = np.arange(len(days_sorted))  # day index
fig = go.Figure(data=[go.Surface(
    x=x_common,   # log return values
    y=Y,          # day index
    z=Z,          # densities
    colorscale='Turbo'
)])
fig.update_layout(
    # title='Daily Kernel Density Estimates of Log Returns',
    scene=dict(
        xaxis_title='Log return',
        yaxis_title='Day',
        zaxis_title='Density',
        yaxis=dict(
            tickmode='array',
            tickvals=np.linspace(0, len(days_sorted)-1, 10, dtype=int),
            ticktext=[str(days_sorted[i]) for i in np.linspace(0, len(days_sorted)-1, 10, dtype=int)]
        )
    ),
    height=700,
    scene_camera=dict(
        eye=dict(x=0.7, y=2.0, z=0.8)  # change these values to rotate
    ),
    title={
        'text': "Daily Kernel Density Estimates of Log Returns",
        'x': 0.5,  # Center horizontally
        'xanchor': 'center'
    }
)
fig.write_image("img/kde_surface.png", width=1200, height=800, scale=2)