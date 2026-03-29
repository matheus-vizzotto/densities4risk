import numpy as np

def inner_product(f, g, du):
    """
    Compute the L^2 inner product ⟨f, g⟩ = ∫ f(u) g(u) du
    in its discrete approximation.

    Parameters
    ----------
    f : ndarray (m, k1)
        First function(s) evaluated on a grid.
    g : ndarray (m, k2)
        Second function(s) evaluated on the same grid.
    du : float
        Grid spacing.

    Returns
    -------
    ndarray (k1, k2)
        Matrix of inner products.
    """
    return np.dot(f.T, g) * du


def L2norm(f, du):
    """
    Compute the L^2 norm of a function vector f.

    Parameters
    ----------
    f : ndarray (m,)
        Function values on a grid.
    du : float
        Grid spacing.

    Returns
    -------
    float
        The L^2 norm.
    """
    return np.sqrt(inner_product(f, f, du))



############################
######### Wavelets #########
############################

def wavedec_sizes(signal_length, wavelet_name, level):
    """
    Compute expected pywt.wavedec sizes:
      returns sizes = [len(A_N), len(D_N), ..., len(D_1)]
      and total = sum(sizes)
    """
    w = pywt.Wavelet(wavelet_name)
    F = w.dec_len
    L = signal_length

    lengths = []
    Lj = L
    for j in range(level):
        Lj = (Lj + F - 1) // 2
        lengths.append(Lj)

    sizes = [lengths[-1]] + lengths[::-1]

    total = sum(sizes)
    return sizes, total

# split wavelet vector into coeff lists
def vec2coeffs(vec, sizes):
    coeffs = []
    pos = 0
    for s in sizes:
        coeffs.append(vec[pos:pos+s])
        pos += s
    return coeffs