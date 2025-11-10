import numpy as np

def inner_product(f, g, du):
    """
    f, g = vectors
    """
    ip = (f.T @ g)*du
    return ip

def L2norm(f,du):
    norm = np.sqrt(inner_product(f, f, du))
    return norm

def sampleCols(A):
    """
    Randomly shuffle columns of a matrix.
    
    Parameters:
    A : numpy.ndarray
        Input matrix (2D array)
    
    Returns:
    numpy.ndarray : Matrix with columns randomly shuffled
    """
    # Get number of columns
    n_cols = A.shape[1]
    
    # Generate random permutation of column indices
    idx = np.random.permutation(n_cols)
    
    # Create new matrix with shuffled columns
    M = A[:, idx].copy()
    
    return M

def super_fun(
        Y,
        lag_max,
        B,
        alpha,
        du,
        p,
        m,
        u,
        select_ncomp,
        dimension
        ) -> dict:
    n = N = Y.shape[1]
    Ybar = np.mean(Y, axis=1, keepdims=True)
    Ydev = Y - Ybar
    core = inner_product(Ydev, Ydev, du=0.05)
    Kstar_core0 = core[:(n-p),:(n-p)]
    Kstar_core = np.zeros((n-p, n-p, p))
    for k in range(p):
        Kstar_core[:, :, k] = core[k : n-(p-k)-1, k : n-(p-k)-1]
    Kstar_sum = np.zeros([(n-p),(n-p)]) 
    for k in range(p):
        Kstar_sum = Kstar_sum + Kstar_core[:,:,k]
    Kstar = ((n-p)**(-2)) * (Kstar_sum.T @ Kstar_core0)
    eigenvalues, eigenvectors = np.linalg.eig(Kstar)
    idx = eigenvalues.argsort()[::-1] 
    thetahat = eigenvalues[idx]
    gammahat = eigenvectors[:,idx]
    tol = 10**(-4)
    for j in range(np.min([len(thetahat),10])):
        if (abs(thetahat[j].imag)>tol):
            print("Complex eigenvalue found.")
    thetahat = thetahat.real
    gammahat = gammahat.real
    thetahat_old = thetahat
    gammahat_old = gammahat

    if select_ncomp:
        bs_pvalues = np.full(lag_max, 0)
        for d0 in range(lag_max-1): 
            thetahatH0 = thetahat_old[d0+1]
            thetahat   = thetahat_old[d0]
            gammahat   = gammahat_old[:,:(d0+1)]
            psihat_root = Ydev[:,:(n-p)] @ gammahat_old[:,:(d0+1)]
            psihat = np.zeros([m, d0+1])
            for i in range(d0+1):
                print(i)
                psihat[ : , i ] = psihat_root[:,i]/L2norm(psihat_root[:,i], du=du)
            etahat = inner_product(psihat, Ydev, du=du)
            Yhat = Ybar + (psihat @ etahat) 
            Yhat_fix = Yhat
            if (Yhat_fix[Yhat_fix<0]).size != 0:
                print(Yhat_fix)
                print("Negative values were found; substituting them with 0.")
                Yhat_fix[np.where(Yhat_fix < 0)] = 0
                for t in np.range(N):
                    Yhat_fix[:,t] = Yhat_fix[:,t]/(sum(Yhat_fix[:,t])*du) 
            epsilonhat = Y - Yhat_fix

            bs_thetahat = np.full(B,0)
            bs = {}
            for i in range(B):
                bs["epsilon"] = sampleCols(epsilonhat)
                bs["Y"]       = Yhat_fix + bs["epsilon"]
                bs["Ybar"]    = np.mean(bs["Y"], axis=1, keepdims=True)
                bs["Ydev"]    = bs["Y"] - bs["Ybar"]
                bs["core"]    = inner_product(bs["Ydev"], bs["Ydev"], du = du)
                bs["Kstar_core0"] = bs["core"][:(n-p),:(n-p)]
                bs["Kstar_core"]  = np.zeros((n-p, n-p, p))
                for k in range(p):
                    bs["Kstar_core"][:,:,k] = bs["core"][(k+1):(n-(p-k)),(k+1):(n-(p-k))]
                bs["Kstar_sum"] = np.zeros([n-p,n-p])
                for k in range(p):
                    bs["Kstar_sum"] = bs["Kstar_sum"] + bs["Kstar_core"][:,:,k]
                bs["Kstar"] = (n-p)**(-2) * bs["Kstar_sum"] @ bs["Kstar_core0"]
                eigenvalues, _ = np.linalg.eig(bs["Kstar"])
                idx = eigenvalues.argsort()[::-1] 
                bs["thetahat"] = eigenvalues[idx]
                bs_thetahat[i] = bs["thetahat"][d0+1]
            bs_pvalues[d0] = (bs_thetahat >= thetahatH0).sum()/B
        d0 = np.where(bs_pvalues < alpha)[0][0] + 1
        if np.isnan(d0):
            d0 = 1
    else:
        d0 = dimension
        
    thetahat = (thetahat_old[:d0+1]).real
    gammahat = (gammahat_old[:,:d0+1]).real
    psihat_root = Ydev[:,:(n-p)] @ gammahat
    psihat = np.zeros([m,d0+1])
    for i in range((d0+1)):
        psihat[:,i] = psihat_root[:,i]/L2norm(psihat_root[:,i],du)
    etahat = inner_product(psihat, Ydev, du)
    Yhat = Ybar + psihat @ etahat
    Yhat_fix = Yhat
    if (Yhat_fix[Yhat_fix<0]).size != 0:
        print(Yhat_fix)
        print("Negative values were found; substituting them with 0.")
        Yhat_fix[np.where(Yhat_fix < 0)] = 0
        for t in np.range(N):
            Yhat_fix[:,t] = Yhat_fix[:,t]/(sum(Yhat_fix[:,t])*du)
    epsilonhat = Y - Yhat_fix

    result = {
        "Y": Y,
        "Ybar": Ybar,
        "thetahat": thetahat,
        "gammahat": gammahat,
        "psihat": psihat,
        "etahat": etahat,
        "Yhat": Yhat,
        "Yhat_fix": Yhat_fix,
        "epsilonhat": epsilonhat,
        "u": u,
        "d0": d0
        }    
    if select_ncomp:
        result["bs_pvalues"] = bs_pvalues

    return result