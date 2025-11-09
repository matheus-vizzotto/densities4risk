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

Ybar = rowMeans(Y)
Ydev = Y - Ybar

core = inner.product(Ydev,Ydev, du)
core
Kstar.core0 = core[1:(n-p),1:(n-p)]
Kstar.core0


Kstar.core = array(0,c(n-p,n-p,p))
Kstar.core