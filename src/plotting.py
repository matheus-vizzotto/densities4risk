import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
from typing import List
from plotly.subplots import make_subplots
from statsmodels.tsa.stattools import acf, pacf

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

def set_plotting_configs():
    plt.rcParams.update({

        # ---- Figure layout ----
        "figure.figsize": (12, 5),       # Nice default for time-series
        "figure.dpi": 110,
        "axes.spines.top": False,
        "axes.spines.right": False,

        # ---- Background colors ----
        "axes.facecolor": "#f7f7f7",
        "figure.facecolor": "#ffffff",
        "savefig.facecolor": "#ffffff",

        # ---- Grids ----
        "axes.grid": True,
        "grid.color": "0.85",
        "grid.linestyle": "--",
        "grid.linewidth": 0.7,

        # ---- Fonts ----
        "font.size": 12,
        "axes.titlesize": 15,
        "axes.labelsize": 13,

        # ---- Lines ----
        "lines.linewidth": 1,

        # ---- Ticks ----
        "xtick.major.size": 0,
        "ytick.major.size": 0,

        # ---- Legend ----
        "legend.frameon": False,
    })

    # Optional: nicer color cycle (using seaborn palette)

    import seaborn as sns
    colors = sns.color_palette("deep")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=colors)


def px_select_menu(fig):
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=list([
                    dict(
                        args=["visible", "legendonly"],
                        label="Deselect All",
                        method="restyle"
                    ),
                    dict(
                        args=["visible", True],
                        label="Select All",
                        method="restyle"
                    )
                ]),
                pad={"r": 10, "t": 10},
                showactive=False,
                x=1,
                xanchor="right",
                y=1.1,
                yanchor="top"
            ),
        ]
    )

import re

def style_modes_of_variation(
        fig, 
        alphas: List[float], 
        mu_transparency: float=0.2,
        mu_black:bool=True
        ):
    """
    Styles a Plotly figure to visualize FPCA modes of variation with mathematical labels.
    
    Positive perturbations are colored green and negative perturbations are red. 
    The 'mean' (alpha=0) is styled as a semi-transparent thick black line to act as a 
    baseline. Line styles (dotted to solid) are assigned based on the absolute 
    intensity of alpha to visually represent the magnitude of the shift.

    Parameters:
    ----------
    fig : plotly.graph_objects.Figure
        The figure object generated by plotly express (long format).
    alphas : list of float
        The list of alpha intensities used in the perturbation (e.g., [0, 0.5, 2, 10]).
    mu_transparency : float, default 0.2
        The alpha channel (0 to 1) for the black mean line. Lower is more transparent.

    Returns:
    -------
    plotly.graph_objects.Figure
        The styled figure with updated legend names, colors, and hover templates.
    """
    # 1. Line styles from least to most continuous (matches increasing intensity)
    dash_styles = ["dot", "dashdot", "dash", "longdash"]
    
    # 2. Extract unique absolute intensities (excluding 0)
    sorted_abs_alphas = sorted(list(set([abs(float(a)) for a in alphas if float(a) != 0])))
    
    # 3. Map intensities to styles (higher alpha = more continuous line)
    style_map = {val: (dash_styles[i] if i < len(dash_styles) else "solid") 
                 for i, val in enumerate(sorted_abs_alphas)}

    def style_traces(trace):
        # --- MEAN CONFIGURATION & LABEL ---
        if trace.name in ["alpha_0", "alpha_0.0", "mean"]:
            if mu_black:
                mu_color = f"rgba(0, 0, 0, {mu_transparency})"
            else:
                mu_color = f"rgba(255, 255, 255, {mu_transparency})"
            trace.line.update(color=mu_color, width=5, dash="solid")
            trace.name = "alpha=0"
            trace.hovertemplate = "<b>alpha=0 (Mean)</b><br>x: %{x}<br>y: %{y}<extra></extra>"
            return

        # --- EXTRACT VALUE & DIRECTION ---
        match = re.search(r"alpha_([\d\.]+)", trace.name)
        if match:
            val_float = float(match.group(1))
            
            if "minus" in trace.name:
                final_val = -val_float
                trace.line.color = "red"
            else:
                final_val = val_float
                trace.line.color = "green"
            
            # Formatting: alpha=10 instead of alpha=10.0
            label_val = int(final_val) if final_val == int(final_val) else final_val
            trace.name = f"alpha={label_val}"

            # --- ASSIGN LINE STYLE ---
            trace.line.dash = style_map.get(abs(val_float), "solid")
            
            # --- UPDATE HOVER LABEL ---
            trace.hovertemplate = f"<b>alpha={label_val}</b><br>x: %{{x}}<br>y: %{{y}}<extra></extra>"

    fig.for_each_trace(style_traces)
    
    # Ensure legend is organized and hover is unified for easy comparison
    fig.update_layout(
        legend={'traceorder': 'grouped'},
        hovermode="x unified"
    )
    
    return fig




import numpy as np
from statsmodels.tsa.stattools import acf, pacf
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def plot_acf_pacf(
    series,
    nlags=10,
    alpha=0.05,
    title=None,
):
    """
    Side-by-side ACF and PACF plots.
    """

    name = getattr(series, "name", "Series")

    y = np.asarray(series).astype(float)
    y = y[~np.isnan(y)]

    # Subplots with simple titles
    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("ACF", "PACF"),
        horizontal_spacing=0.12
    )

    # -------- ACF --------
    acf_array = acf(y, alpha=alpha, nlags=nlags)
    acf_vals = acf_array[0]
    acf_lower = acf_array[1][:, 0] - acf_vals
    acf_upper = acf_array[1][:, 1] - acf_vals

    lags = np.arange(len(acf_vals))

    for x in lags:
        fig.add_scatter(
            x=(x, x),
            y=(0, acf_vals[x]),
            mode="lines",
            line_color="#3f3f3f",
            row=1,
            col=1
        )

    fig.add_scatter(
        x=lags,
        y=acf_vals,
        mode="markers",
        marker_color="#1f77b4",
        marker_size=9,
        row=1,
        col=1
    )

    fig.add_scatter(
        x=lags,
        y=acf_upper,
        mode="lines",
        line_color="rgba(255,255,255,0)",
        row=1,
        col=1
    )

    fig.add_scatter(
        x=lags,
        y=acf_lower,
        mode="lines",
        fill="tonexty",
        fillcolor="rgba(32,146,230,0.25)",
        line_color="rgba(255,255,255,0)",
        row=1,
        col=1
    )

    # -------- PACF --------
    pacf_array = pacf(y, alpha=alpha, nlags=nlags)
    pacf_vals = pacf_array[0]
    pacf_lower = pacf_array[1][:, 0] - pacf_vals
    pacf_upper = pacf_array[1][:, 1] - pacf_vals

    lags = np.arange(len(pacf_vals))

    for x in lags:
        fig.add_scatter(
            x=(x, x),
            y=(0, pacf_vals[x]),
            mode="lines",
            line_color="#3f3f3f",
            row=1,
            col=2
        )

    fig.add_scatter(
        x=lags,
        y=pacf_vals,
        mode="markers",
        marker_color="#1f77b4",
        marker_size=9,
        row=1,
        col=2
    )

    fig.add_scatter(
        x=lags,
        y=pacf_upper,
        mode="lines",
        line_color="rgba(255,255,255,0)",
        row=1,
        col=2
    )

    fig.add_scatter(
        x=lags,
        y=pacf_lower,
        mode="lines",
        fill="tonexty",
        fillcolor="rgba(32,146,230,0.25)",
        line_color="rgba(255,255,255,0)",
        row=1,
        col=2
    )

    # Axis styling
    fig.update_yaxes(zerolinecolor="#000000")
    fig.update_traces(showlegend=False)

    if title is None:
        title = f"ACF & PACF of {name}"

    fig.update_layout(
        title=title,
        title_x=0.5,
        template="plotly_white",
        height=450
    )

    return fig

def plot_ccfs(
    series1,
    series2,
    nlags=20,
    alpha=0.05,
    title=None,
):
    """
    Cross-correlation plot between two series with symmetric lags.
    """

    name1 = getattr(series1, "name", "Series 1")
    name2 = getattr(series2, "name", "Series 2")

    y1 = np.asarray(series1).astype(float)
    y2 = np.asarray(series2).astype(float)

    mask = ~np.isnan(y1) & ~np.isnan(y2)
    y1, y2 = y1[mask], y2[mask]

    # Standardize (important for CCF)
    y1 = (y1 - y1.mean()) / y1.std()
    y2 = (y2 - y2.mean()) / y2.std()

    n = len(y1)

    # Full cross-correlation
    corr_full = np.correlate(y1, y2, mode="full") / n
    mid = len(corr_full) // 2

    corr = corr_full[mid - nlags : mid + nlags + 1]
    lags = np.arange(-nlags, nlags + 1)

    # Asymptotic confidence bands
    z = abs(np.quantile(np.random.standard_normal(200000), 1 - alpha / 2))
    ci = z / np.sqrt(n)

    upper = np.full_like(corr, ci)
    lower = -upper

    # -------- Plot --------
    fig = go.Figure()

    for i, x in enumerate(lags):
        fig.add_scatter(x=(x, x), y=(0, corr[i]),
                        mode="lines", line_color="#3f3f3f")

    fig.add_scatter(x=lags, y=corr, mode="markers",
                    marker_color="#1f77b4", marker_size=10)

    fig.add_scatter(x=lags, y=upper, mode="lines",
                    line_color="rgba(255,255,255,0)")

    fig.add_scatter(x=lags, y=lower, mode="lines",
                    fill="tonexty", fillcolor="rgba(32,146,230,0.3)",
                    line_color="rgba(255,255,255,0)")

    fig.update_traces(showlegend=False)
    fig.update_yaxes(zerolinecolor="#000000")

    if title is None:
        title = f"Cross-Correlation of {name1} × {name2}"

    fig.update_layout(title=title)
    fig.show()

import numpy as np
import pandas as pd
import plotly.graph_objects as go


def plot_density_timeseries_3d(
    df: pd.DataFrame,
    x_range=None,
    colorscale="Viridis",
    camera=None,
    theme="plotly_white",
    percentile_cut: List[float] = None,
    title = None
):
    """
    3D surface plot of density time series.

    Parameters
    ----------
    df : pd.DataFrame
        First column = index/grid (x)
        Remaining columns = dates (y)
        Values = pdf (z)

    x_range : tuple or None, optional
        (xmin, xmax) to restrict x-axis domain.

    colorscale : str
        Plotly colorscale name (e.g. 'Viridis', 'Cividis', 'Plasma', 'Turbo').

    camera : dict or None
        Plotly camera dictionary, e.g.
        dict(eye=dict(x=-1.5, y=1.5, z=1.0))

    theme : str
        Plotly template (e.g. 'plotly_white', 'plotly_dark', 'ggplot2').
    """

    x = df.index

    if x_range is not None:
        xmin, xmax = x_range
        mask = (x >= xmin) & (x <= xmax)
        df = df.loc[mask]
        x = x[mask]

    y = pd.to_datetime(df.columns[1:])
    z = df.iloc[:, 1:].T.values

    if percentile_cut is None:
        cmin = np.percentile(z,0)
        cmax = np.percentile(z,100)
    else:
        cmin = np.percentile(z,percentile_cut[0])
        cmax = np.percentile(z,percentile_cut[1])

    fig = go.Figure(
        data=[
            go.Surface(
                x=x,
                y=y,
                z=z,
                colorscale=colorscale,
                cmin=cmin,
                cmax=cmax
            )
        ]
    )

    fig.update_layout(
        title=title,
        template=theme,
        scene=dict(
            xaxis_title="Index",
            yaxis_title="Date",
            zaxis_title="PDF",
            xaxis=dict(autorange="reversed"),
        ),
        height=700,
    )

    if camera is not None:
        fig.update_layout(scene_camera=camera)

    return fig


