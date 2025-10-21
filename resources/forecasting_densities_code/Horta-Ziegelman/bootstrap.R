####################################################################
## SEQUENTIAL BOOTSTRAP for selecting the number of components d ##
####################################################################

################################################################################
# Defining a function sampleCols, which samples the columns of any given matrix
# Randomly sampling the column index
################################################################################

sampleCols = function(A)
{
    idx = sample(1:ncol(A))
    M = matrix(0, nrow = nrow(A), ncol = ncol(A))
    for(i in 1:length(idx))
    {
        M[,i] = A[,idx[i]]
    }
    return(M)
}

# This is a sequential bootstrap which stores the p-value associated to the (d0+1)-th
# eigenvalue, with d0 running from 1 through 2.
# The

# Vector which stores the p-values
bs.pvalues = numeric()

# Building the estimated functions Yhat with dimension d0, d0 in {1,2,...,lag_max}.
# See the main code for an explanation of this section.

for (d0 in 1:lag_max)
{
    thetahatH0 = Re(thetahat.old[d0+1])
    thetahat = Re(thetahat.old[1:d0])
    gammahat = Re(gammahat.old[,1:d0])
    psihat.root = Ydev[,1:(n-p)] %*% gammahat
    psihat = matrix(0, nrow = m, ncol = d0)
    for(i in 1:d0) psihat[,i] = psihat.root[,i]/L2norm(psihat.root[,i]); rm(i)

    etahat = inner.product(psihat, Ydev)
    Yhat = Ybar + psihat %*% etahat
    epsilonhat = Y - Yhat

    ###############
    ## Bootstrap ##
    ###############
    # Number of bootstrap loops
    B = 100

    # Creating the vector bs.thetahat, which stores the eigenvalues obtained
    # from each bootstrap loop. length(bs.thetahat) = B
    bs.thetahat = numeric()

    bs = list()
    for(i in 1:B)
    {
        # This section of the code repeats c-main.r, however inside the list object bs

        # Building the bootstrap 'observed' function Y
        bs$epsilon = sampleCols(epsilonhat)
        bs$Y = Yhat + bs$epsilon

        bs$Ybar = rowMeans(bs$Y)
        bs$Ydev = bs$Y - bs$Ybar

        bs$core = inner.product(bs$Ydev,bs$Ydev)
        bs$Kstar.core0 = bs$core[1:(n-p),1:(n-p)]
        bs$Kstar.core = array(0,c(n-p,n-p,p))
        for(k in 1:p) bs$Kstar.core[,,k] = bs$core[(k+1):(n-(p-k)),(k+1):(n-(p-k))]
        bs$Kstar.sum = matrix(0,nrow=n-p,ncol=n-p)
        for (k in 1:p) bs$Kstar.sum = bs$Kstar.sum + bs$Kstar.core[,,k]
        bs$Kstar = (n-p)^(-2) * bs$Kstar.sum %*% bs$Kstar.core0
        bs$thetahat = eigen(bs$Kstar)$values
        bs$thetahat = sort(bs$thetahat,decreasing=TRUE)
        # Storing the (d0+1)-th eigenvalue obtained in the i-th bootstrap loop
        bs.thetahat[i] = Re(bs$thetahat[d0+1])
        print(paste('d0 = ', d0, '; Loop = ', i, sep = ""))
    }
    # Storing the p-values
    bs.pvalues[d0] = sum(bs.thetahat >= thetahatH0)/B
}
rm(bs, bs.thetahat, thetahatH0)

# BOOTSTRAP p-values
# 1. h_scale = 1
# > bs.pvalues # p = 1, B = 1000
# [1] 0.00 0.24
# > bs.pvalues # p = 2, B = 1000
# [1] 0.003 0.431
# > bs.pvalues # p = 3, B = 1000
# [1] 0.014 0.730
# > bs.pvalues # p = 4, B = 1000
# [1] 0.037 0.936
# > bs.pvalues # p = 5, B = 1000
# [1] 0.066 0.561
# > bs.pvalues # p = 6, B = 10000
# [1] 0.0362 0.1674
# 
# 2. h_scale = .5
# > bs.pvalues # p = 5, B = 10000
# [1] 0.1308 0.5272
# > bs.pvalues # p = 6, B = 10000
# [1] 0.0130 0.3435
# 
# 3. h_scale = 2
# > bs.pvalues # p = 5, B = 10000
# [1] 0.0102 0.3778
# > bs.pvalues # p = 6, B = 10000
# [1] 0.0269 0.1619
