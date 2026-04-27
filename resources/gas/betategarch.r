install.packages("betategarch")
library("betategarch")


##simulate 1000 observations from model with default parameter values:
set.seed(123)
y <- tegarchSim(1000)
##estimate and store as'mymod':
mymod <- tegarch(y)
##print estimates and standard errors:
print(mymod)
##graph of fitted volatility (conditional standard deviation):
plot(fitted(mymod))
##plot forecasts of volatility 1-step ahead up to 10-steps ahead:
plot(predict(mymod, n.ahead=10))


##generate 1000 random numbers from the skewed t:
eps <- rST(500, df=5) #symmetric t
eps <- rST(500, df=5, skew=0.8) #skewed to the left
eps <- rST(500, df=5, skew=2) #skewed to the right
##compare empirical mean with analytical:
mean(eps)
STmean(5, skew=2)
##compare empirical variance with analytical:
var(eps)
STvar(5, skew=2)


##1-component specification: simulate series with 500 observations:
set.seed(123)
y <- tegarchSim(500, omega=0.01, phi1=0.9, kappa1=0.1, kappastar=0.05,
                df=10, skew=0.8)
##simulate the same series, but with more output (volatility, log-volatility or
##lambda, lambdadagger, u and epsilon)
##plot the simulated values:
plot(y)
##2-component specification: simulate series with 500 observations:
set.seed(123)
y <- tegarchSim(10, omega=0.01, phi1=0.95, phi2=0.9, kappa1=0.01, kappa2=0.05,
                kappastar=0.03, df=10, skew=0.8, verbose=TRUE)

yt     <- coredata(y$y)
sigma  <- coredata(y$sigma)

dens_obs <- (1 / sigma) * dST(yt / sigma, df = 10, skew = 0.8)

grid <- seq(min(yt) - 2, max(yt) + 2, length = 500)

dens_list <- lapply(seq_along(yt), function(t) {
  (1 / sigma[t]) * dST(grid / sigma[t], df = 10, skew = 0.8)
})

dens_mat <- sapply(seq_along(yt), function(t) {
  (1 / sigma[t]) * dST(grid / sigma[t], df = 10, skew = 1)
})

matplot(grid, dens_mat, type = "l", lty = 1,
        main = "Sample of conditional densities",
        xlab = "y", ylab = "Density")
     
grid
