import pandas as pd
import numpy as np
from typing import Tuple, Iterable


def _mode_of_variation(
    mu: pd.Series,
    psi: pd.Series,
    theta: float,
    alpha: float
    ) -> Tuple[pd.Series, pd.Series]:
    """
    Construct plus/minus FPCA modes in L2 space.

    Parameters
    ----------
    mu : pd.Series
        Mean function (in transformed L2 space).

    psi : pd.Series
        Eigenfunction corresponding to a given FPCA component.

    theta : float
        Eigenvalue associated with `psi`.

    alpha : float
        Scaling parameter controlling excursion size
        (typically in multiples of sqrt(theta)).

    Returns
    -------
    (mode_plus, mode_minus) : Tuple[pd.Series, pd.Series]
        Positive and negative perturbations around the mean:
            mu ± alpha * sqrt(theta) * psi
    """

    perturbation = alpha * np.sqrt(theta) * psi

    return mu + perturbation, mu - perturbation


def modes_of_variation(
    df: pd.DataFrame,
    psihat: pd.DataFrame,
    thetahat: Iterable[float],
    alphas: Iterable[float]
    ) -> pd.DataFrame:
    """
    Compute FPCA modes of variation in L2 space.

    This constructs synthetic curves of the form

        mu ± alpha * sqrt(theta_k) * psi_k

    for each principal component k and each scaling factor alpha.

    Parameters
    ----------
    df : pd.DataFrame
        Original transformed functions (columns = time, rows = grid).

    psihat : pd.DataFrame
        Estimated eigenfunctions (rows = grid, columns = components).

    thetahat : iterable of float
        Corresponding eigenvalues.

    alphas : iterable of float
        Multipliers controlling mode magnitude (e.g. [0.5, 1, 2]).

    Returns
    -------
    modes : Long DataFrame with [principal_component, alpha_direction, value]

    Example
    --------
    modes = modes_of_variation(
        df_lqds,
        KdFPC_model.psihat,
        KdFPC_model.thetahat,
        alphas=[0.0, 0.5, 1.0, 1.5, 3.0]
        )

    Notes
    -----
    These modes live in L2 (transform space). To interpret them as
    densities, they must be mapped back via the inverse mLQD transform.
    The inverse should use the mean c=F(0) for all.
    """

    # Mean function in L2
    mu = df.mean(axis=1)

    # List to collect all data rows
    all_records = []

    for alpha in alphas:
        # Modes of variation for each Principal Component
        for k in range(psihat.shape[1]):
            psi_k = psihat[:, k]
            theta_k = thetahat[k]

            if alpha == 0:
                mean, _ = _mode_of_variation(mu, psi_k, theta_k, alpha)
                for x, val in mean.items():
                    all_records.append({
                        'index': x,
                        'pc': f'PC{k+1}',
                        'alpha': f'alpha_{alpha}',
                        'value': val
                    })
            else:
                plus, minus = _mode_of_variation(mu, psi_k, theta_k, alpha)
                
                # Add "plus" variation records
                for x, val in plus.items():
                    all_records.append({
                        'index': x,
                        'pc': f'PC{k+1}',
                        'alpha': f'alpha_{alpha}_plus',
                        'value': val
                    })
                
                # Add "minus" variation records
                for x, val in minus.items():
                    all_records.append({
                        'index': x,
                        'pc': f'PC{k+1}',
                        'alpha': f'alpha_{alpha}_minus',
                        'value': val
                    })

    # Create the DataFrame
    df_long = pd.DataFrame(all_records)
    
    return df_long