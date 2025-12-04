import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go

# TODO: plot_3d_ts traz domínio da função invertido

def plot_2d_fts(
        df : pd.DataFrame, 
        x: np.array, 
        title_="Functional time series", 
        alpha_=0.1,
        output_path=None):
    """
    df (pd.Dataframe): Dataframe with only function values where columns are functional objects and rows are grid points.
    x  (np.array): numpy array of x axis values

    ex: plot_2d_fts(df_densities, x=df_densities.index.values, title_ = "Densities over time")
    """
    cols = df.columns

    plt.figure(figsize=(15,5))
    for col in cols:
        plt.plot(x, df.loc[:, col], c="gray", alpha=alpha_)
    pointwise_mean = df.mean(axis=1)
    plt.plot(x, pointwise_mean, c="red", label="Point-wise mean")

    plt.title(title_)
    plt.legend()

    if not output_path:
        plt.show()
    elif output_path:
        plt.savefig(output_path) # png

def plot_3d_fts(
        df: pd.DataFrame, 
        f_domain: np.array, 
        time_index: np.array, 
        title_="Functional Time Series", 
        axis_titles=["x","y","z"],
        color_scale = "Turbo",
        output_angle = [-0.7, -2.0, 0.8],
        output_path=None):
    """
    ex: plot_3d_fts(df_lqds.iloc[1:-1,:], t[1:-1], df_lqds.columns, title_="LQDensities", axis_titles=["u", "Day", "LQDensity"])
    """
    domain = f_domain
    T = time_index
    image = [df.iloc[:,col] for col in range(df.shape[1])]
    # image = df.to_numpy()

    fig = go.Figure(data=[go.Surface(
        x=domain, 
        y=T,
        z=image,  
        colorscale=color_scale
    )])

    fig.update_layout(
        # title='Daily Kernel Density Estimates of Log Returns',
        scene=dict(
            xaxis_title=axis_titles[0],
            yaxis_title=axis_titles[1],
            zaxis_title=axis_titles[2],
        ),
        height=700,
        scene_camera=dict(
            eye=dict(x=output_angle[0], y=output_angle[1], z=output_angle[2])  # change these values to rotate
        ),
        title={
            'text': title_,
            'x': 0.5,  # Center horizontally
            'xanchor': 'center'
        }
    )

    if not output_path:
        fig.show()
    elif output_path:
        fig.show()
        fig.write_image(output_path, width=1200, height=800, scale=2) # png