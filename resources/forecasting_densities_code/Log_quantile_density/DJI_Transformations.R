# Demonstrate modified LQD transformation for DJI cross-sectional return densities

setwd('/Users/alexanderpetersen/Dropbox/Piotr_DensityTimeSeries/Code_Data')

load('den_v2.rdata')
source('dens2lqdNew.R')
library(fdapace) # has trapezoidal integration used in the transformation functions, also FPCA

data = DJI_ret
kernel = "gaussian"

    # Sample size
    N = nrow(data)

    if (!exists('h_scale')) h_scale = 1
    if(kernel == "gaussian")
    {
        h.hat_5m = sapply(1:N, function(t) 1.06*sd(data[t,])*(length(data[t,])^(-(1/5))))
    }
    if(kernel == "epanechnikov")
    {
        h.hat_5m = sapply(1:N, function(t) 2.34*sd(data[t,])*(length(data[t,])^(-(1/5))))  
    }
    h.hat_5m = h_scale * h.hat_5m

    # 2. Discretization
    # Evaluation points
    u = seq(from = min(data), to = max(data), length = m)
    
    # Interval length
    du = u[2] - u[1]
    
    # Creating an (m x n) matrix which represents the observed densities. Y[j,t] is the density at date t evaluated at u[j]
    if(kernel == "gaussian")
    {
        Y = sapply(1:N, function(t) density(data[t,], bw = h.hat_5m[t], kernel = 'gaussian', from = min(data), to = max(data), n = m)$y)
    }
    if(kernel == "epanechnikov")
    {
        Y = sapply(1:N, function(t) density(data[t,], bw = h.hat_5m[t], kernel = 'epanechnikov', from = min(data), to = max(data), n = m)$y)        
    }

    # correcting to ensure integral Y_t du = 1
    for(t in 1:N)
    {
        Y[,t] = Y[,t]/(sum(Y[,t])*du)
    }

    # Renormalize Densities to have integral 1
    n = ncol(Y)
    N = length(u)

    any(apply(Y, 2, function(z) any(z < 0))) # make sure no density estimates are negative
    # FALSE

    dens = sapply(1:n, function(i){
      Y[,i]/trapzRcpp(X = u, Y = Y[,i])
    })

    # Plot densities

    savepdf("original_density", width = 12, height = 10, toplines = 0.8)
    matplot(u, dens, type = 'l', main = 'Original Density Estimates')
    dev.off()

    # Try Forward transformation
    M = 3000 # number of gridpoints for LQD functions - chosen large here so that 0 isn't too close to the boundary of all supports
    lqd = matrix(0, nrow = M, ncol = n)
    c = rep(0, 1, n)
    t = seq(0, 1, length.out = M)
    t0 = u[which.min(abs(u))] # closest value to 0

    for(i in 1:n)
    {
        tmp = dens2lqd(dens = dens[,i], dSup = u, lqdSup = t, t0 = t0, verbose = FALSE)
        lqd[,i] = tmp$lqd
        c[i] = tmp$c
    }

    ## Think about logit transformation or log ratio

    # Plot LQD Functions

    savepdf("LQD", width = 22, height = 10, toplines = 0.8)
    par(mfrow = c(1,2))
    matplot(t, lqd, type = 'l', main = 'LQD functions') # Shows problems near boundaries
    matplot(t, lqd, type = 'l', ylim = c(min(lqd[lqd < Inf]), 10), main = 'LQD functions')
    dev.off()

    # Now Try Backward Transform

    source('lqd2densNew.R')

    # The following results are useless - they demonstrate the instability of the inverse transformation when the LQD is large at the boundary
    res1 = lapply(1:n, function(i) lqd2dens(lqd = lqd[,i], lqdSup = t, t0 = t0, c = c[i], useSplines = FALSE, verbose = FALSE)) # minimal boundary cuts
    dens1 = sapply(res1, function(r){
        approx(x = r$dSup, y = r$dens, xout = u, yleft = 0, yright = 0)$y})
    matplot(u, dens1, type = 'l', main = 'Density Estimates after Transformation - No Boundary Cutting')

    # Try cutting off boundary points with large LQD values (effectively setting the density to zero rather than trying to compute it numerically)
    cut = list()
    res2 = list()
    for(i in 1:n){ 
      cut[[i]] = c(0, 0)
      cut[[i]][1] = sum(lqd[1:15,i] > 7)
      cut[[i]][2] = sum(lqd[(M-14):M, i] > 7)
      res2[[i]] = lqd2dens(lqd = lqd[,i], lqdSup = t, t0 = t0, c = c[i], cut = cut[[i]], useSplines = FALSE, verbose = FALSE)
    }

    dens2 = sapply(res2, function(r){
        approx(x = r$dSup, y = r$dens, xout = u, yleft = 0, yright = 0)$y
    })
    
    savepdf("Density_recon", width = 24, height = 10, toplines = 0.8)
    par(mfrow = c(1, 2))
    matplot(u, dens, type = 'l', main = 'Original Density Estimates')
    matplot(u, dens2, type = 'l', main = 'Density Estimates after Forward/Backward Transform')
    dev.off()

    # Assess loss of mass incurred by boundary cutoff
    totalMass = apply(dens2, 2, function(d) trapzRcpp(X = u, Y = d))
    min(totalMass) # 0.9993
    max(totalMass) # 0.99998

    # Numerical comparison of densities 
    
    L2Diff = sapply(1:n, function(i) sqrt(trapzRcpp(X = u, Y = (dens[,i] - dens2[,i])^2))) # L^2 norm difference
    unifDiff = sapply(1:n, function(i){
      interior = which(dens2[,i] > 0)
      max(abs(dens[interior,i] - dens2[interior,i])) 
    }) # Uniform Metric excluding missing boundary values (due to boundary cutoff)
