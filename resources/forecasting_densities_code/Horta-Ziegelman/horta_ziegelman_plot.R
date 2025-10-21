# Setting up the functions

# Defining the inner product on L^2([a,b])
inner.product = function(f, g, du) drop(crossprod(f, g))*du

# Defining the L2-norm
L2norm = function(f, du) sqrt(inner.product(f, f, du))

# Some functions to calculate the moments of a distribution:
# Moments
m1   = function(f) inner.product(u, f, du)
m2   = function(f) inner.product((u - m1(f))^2, f, du)
m3   = function(f) inner.product((u - m1(f))^3, f, du)
skw  = function(f) m3(f)/(m2(f)^(3/2))
m4   = function(f) inner.product((u - m1(f))^4, f, du)
kurt = function(f) m4(f)/(m2(f)^2) - 3

# Importing the data
setwd("/Volumes/MY PASSPORT/Dropbox/Todos/cs/data")
data = read.table("IBovespa_5_min.csv",header=TRUE,sep=";")
data$Time = as.POSIXct(paste(data$Data, data$Hora), format="%d/%m/%Y %H:%M")
data.list = split(data[, c('Time','Fech')], format(data$Time, format = "%Y%m%d"))

# Sample size
N = length(data.list)

# Creating the list 'returns'; returns[[t]](j) is the j-th return on day t

# 5-minute returns
ret = lapply(1:N,
             function(t) {m = length(data.list[[t]]$Fech); sapply(1:(m-1),
                                                                  function(j) log(data.list[[t]]$Fech[j+1]) - log(data.list[[t]]$Fech[j]))})

# Creating the 'rule of thumb' bandwidths h.hat, to be used to estimate the daily densities. h.hat is a n-vector with h.hat[t] being Silverman's rule of thumb bandwidth for day t
if (!exists('h_scale')) h_scale = 1
h.hat_5m = sapply(1:N,
                  function(t) 1.06*sd(ret[[t]])*(length(ret[[t]])^(-(1/5))))
#function(t) 2.34*sd(ret[[t]])*(length(ret[[t]])^(-(1/5)))) # Epanechnikov rule-of-thumb bw
h.hat_5m = h_scale * h.hat_5m

# Initialization parameters
n = N # Number of daily observations
p = 5   # [See Bathia, Yao, Ziegelmann; page 5)
b = .06
a = -b

# 2. Discretization
# Evaluation points
m = 5001 # number of grid points (this is arbitrary)
u = seq(from = a, to = b, length = m)

# Interval length
du = u[2] - u[1]

# Creating an (m x n) matrix which represents the observed densities. Y[j,t] is the density at date t evaluated at u[j]
Y = sapply(1:N, function(t) density(ret[[t]], bw = h.hat_5m[t], kernel = 'gaussian', from = a, to = b, n = m)$y)

# correcting to ensure integral Y_t du = 1
for(t in 1:N)
{
    Y[,t] = Y[,t]/(sum(Y[,t])*du)
}

source("super_fun.R")
foo_out = super_fun(Y = Y, lag_max = 4, B = 1000, alpha = 0.05, du = du, p = 5, m = 5001, u = u)

# first curve in the sample + "filtered" estimate:
plot(foo_out$u, foo_out$Y[,1], type = 'l', xlim = c(-.005,.005), col = rgb(.7,.7,.7), ann = FALSE)
lines(foo_out$u, foo_out$Yhat[,1], col = rgb(.8,.4,0))

# mean curve + eigenfunctions:
plot(foo_out$u, foo_out$Ybar, type = 'l', xlim = c(-.005,.005), ylim = c(-450,450), ann = FALSE)
lines(foo_out$u, 10 * foo_out$psihat[,1], xlim = c(-.005,.005), col = rgb(0,.4,.8))
lines(foo_out$u, 10 * foo_out$psihat[,2], xlim = c(-.005,.005), col = rgb(.8,.4,0))
# lines(foo_out$u, 10 * foo_out$psihat[,3], xlim = c(-.005,.005), col = rgb(.8,.8,.8), lty = 'dotted') # from the third onwards the eigenfunctions are "ugly"

# log-eigenvalues:
plot(1:10, log(foo_out$thetahat[1:10]), pch = 19, ann = FALSE)

# latent time series (the first 4 of them)
par(mfrow = c(2,2))
for (t in 1:4) plot(etahat[t,], type = 'l', ann = FALSE)
par(mfrow = c(1,1))
