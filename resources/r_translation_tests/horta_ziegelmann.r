v1 <- c(1, 2, 3, 4, 5)
v2 <- c(6, 7, 8, 9, 10)
crossprod(v1, v2)

# Defining the inner product on L^2([a,b])
inner.product = function(f, g, du) {drop(crossprod(f, g))*du}

# Defining the L2-norm
L2norm = function(f, du) {sqrt(inner.product(f, f, du))}

# Create rectangular matrices (3x2 and 3x4)
matrixX <- matrix(c(1, 2,
                   3, 4,
                   5, 6), 
                 nrow = 3, ncol = 2, byrow = TRUE)

Y <- matrix(c(2, 1, 3, 4,
                   1, 0, 2, 1,
                   3, 2, 1, 0,
                   1, 2, 3, 4,
                   4, 3, 2, 1
                   ), 
                 nrow = 5, ncol = 4, byrow = TRUE)

n = N = ncol(Y)
p=5
p=2
du=0.005
m=5001

Ybar = rowMeans(Y)
Ydev = Y - Ybar

core = inner.product(Ydev,Ydev, du)
core
Kstar.core0 = core[1:(n-p),1:(n-p)]
Kstar.core0


Kstar.core = array(0,c(n-p,n-p,p))
Kstar.core

for (k in 1:p) Kstar.core[,,k] = core[(k+1):(n-(p-k)),(k+1):(n-(p-k))]
Kstar.core

# Summing the matrices in 'Kstar.core'
Kstar.sum = matrix(0,nrow=n-p,ncol=n-p)
Kstar.sum
for (k in 1:p) Kstar.sum = Kstar.sum + Kstar.core[,,k]
Kstar.sum


Kstar.sum %*% Kstar.core0

Kstar = (n-p)^(-2) * Kstar.sum %*% Kstar.core0

thetahat = eigen(Kstar)$values
thetahat

gammahat = eigen(Kstar)$vectors
gammahat

tol = 10^(-4)


for (j in 1:10){
  if (abs(Im(thetahat[j]))>tol) print("Complex eigenvalue found.")
}


thetahat = Re(thetahat)
thetahat = sort(thetahat, index.return=TRUE, decreasing=TRUE)
thetahat.index = thetahat$ix
thetahat = thetahat$x


gammahat.temp = matrix(0, nrow=nrow(gammahat), ncol=ncol(gammahat))
for(j in 1:(length(thetahat)))
{
  gammahat.temp[,j] = gammahat[,thetahat.index[j]]
}
gammahat = gammahat.temp



# Storing the original eigenvalues and eigenvectors
thetahat.old = thetahat
gammahat.old = gammahat

lag_max = 3
bs.pvalues = vector("numeric", lag_max)

# if(select_ncomp == "TRUE")
# {
  bs.pvalues = vector("numeric", lag_max)
  
  # Building the estimated functions Yhat with dimension d0, d0 in {1,2,...,lag_max}.
  # See the main code for an explanation of this section.
  
  for(d0 in 1:lag_max)
  {
    thetahatH0 = Re(thetahat.old[d0+1])
    thetahat = Re(thetahat.old[1:d0])
    gammahat = Re(gammahat.old[,1:d0])
    psihat.root = Ydev[,1:(n-p)] %*% gammahat
    psihat = matrix(0, nrow = m, ncol = d0)
    # for(i in 1:d0) psihat[,i] = psihat.root[,i]/L2norm(psihat.root[,i], du = du); rm(i)
  }
# }