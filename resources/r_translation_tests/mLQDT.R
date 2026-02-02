install.packages("fdapace")
install.packages("ftsa")


library("ftsa")
library("fdapace") # -> trapzRcpp 
data(DJI_return)
data <- DJI_return

dens2lqd = function(dens, dSup, lqdSup = seq(0, 1, length.out = length(dSup)), t0 = dSup[1], verbose = TRUE){
  
  # Check density requirements
  if(any(dens < 0)){
    stop('Please correct negative density values.')
  }
  
  if(abs( trapzRcpp(X = dSup, Y = dens) - 1) > 1e-5){
    
    warning('Density does not integrate to 1 with tolerance of 1e-5 - renormalizing now.')
    dens = dens/trapzRcpp(X = dSup, Y = dens)
    
  }
  
  if(any(dens == 0)){
    if(verbose){
      print("There are some zero density values - truncating support grid so all are positive")
    }
    lbd = min(which(dens > 0))
    ubd = max(which(dens > 0))
    dens = dens[lbd:ubd]
    dSup = dSup[lbd:ubd]
    dens = dens/trapzRcpp(X = dSup, Y = dens)
  }
  
  N = length(dSup)
  
  # Check LQD output grid
  if(is.null(lqdSup)){
    
    lqdSup = seq(0, 1, length.out = N)
    
  }else if(!all.equal( range(lqdSup),c(0,1) )){
    
    if(verbose){
      print("Problem with support of the LQD domain's boundaries - resetting to default.")
    }
    lqdSup = seq(0, 1, length.out = N)
    
  }
  
  # Check t0
  if(!(t0 %in% dSup)){
    
    if(verbose){
      print("t0 is not a value in dSup - resetting to closest value")
    }
    t0 = dSup[which.min(abs(dSup - t0))]
    
  }
  
  M = length(lqdSup) 
  c_ind = which(dSup == t0)
  
  # Get CDF and lqd on temporary grid, compute c
  tmp = cumtrapzRcpp(X = dSup, dens)
  c = tmp[c_ind]
  
  indL = duplicated(tmp[1:floor(N/2)])
  indR = duplicated(tmp[(floor(N/2)+1):N], fromLast = TRUE)
  qtemp = tmp[!c(indL, indR)]
  lqd_temp = -log(dens[!c(indL, indR)]);
  
  # Interpolate lqdSup, keeping track of Inf values at boundary, then compute c
  lqd = rep(0, 1, M)
  
  if(any(is.infinite(lqd_temp[c(1, N)]))){
    
    tmpInd = 1:N
    Ind = 1:M
    
    if(lqd_temp[1] == Inf){
      
      lqd[1] = Inf
      tmpInd = tmpInd[-1]
      Ind = Ind[-1]
      
    }
    
    if(lqd_temp[N] == Inf){
      
      lqd[M] = Inf
      tmpInd = tmpInd[-length(tmpInd)]
      Ind = Ind[-length(Ind)]
      
    }
    
    lqd[Ind] = approx(x = qtemp[tmpInd], y = lqd_temp[tmpInd], xout = lqdSup[Ind], rule = 2)$y 
    
  }else{
    
    lqd = approx(x = qtemp, y = lqd_temp, xout = lqdSup, rule = c(2,2))$y 
    
  }
  
  return(list('lqdSup',  lqdSup, 'lqd' = lqd, 'c' = c))
}

kernel = "gaussian"
m = 5001
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

matplot(u, dens, type = 'l', main = 'Original Density Estimates')

i = 100
plot(u, dens[,i], , type="l",lwd=0.5)

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

matplot(t, lqd, type = 'l', main = 'LQD functions') # Shows problems near boundaries
matplot(t, lqd, type = 'l', ylim = c(min(lqd[lqd < Inf]), 10), main = 'LQD functions')

# Backwards transformation
lqd2dens = function(lqd, lqdSup = seq(0, 1, length.out = length(lqd)), dSup, t0 = 0, c = 0, useSplines = TRUE, cut = c(0, 0), verbose = TRUE){
  
  if(!all.equal( range(lqdSup),c(0,1) )){
    
    warning("Problem with support of the LQD domain's boundaries - resetting to default.")
    lqdSup = seq(0, 1, length.out = length(lqd))
    
  }
  
  M = length(lqd)
  r = which(exp(lqd) == Inf)
  
  if(length(r) > 0){
    
    if(any(r < floor(M/2))){
      cut[1] = max(cut[1], max(r[r < floor(M/2)]))
    }
    if(any(r >= floor(M/2))){
      cut[2] = max(cut[2], M - min(r[r >= floor(M/2)]) + 1)
    }    
  }
  
  # Cut boundaries
  lqdSup = lqdSup[(cut[1] + 1):(M - cut[2])]
  lqd = lqd[(cut[1] + 1):(M - cut[2])]
  M = length(lqd) # reset N
  
  if(!(c %in% lqdSup)){
    
    if(c < lqdSup[1] || c > lqdSup[M]){
      
      stop("c is not contained withing range of lqdSup after cutoff")
      
    }
    
    if(verbose){
      
      print("c is not equal to a value in lqdSup - resetting to closest value")
      
    }
    c = lqdSup[which.min(abs(lqdSup - c))]
    
  }
  
  c_ind = which(lqdSup == c)
  
  if( useSplines ){    # Could fit spline if this yields more accurate numerical integration
    
    lqd_sp = splinefun(lqdSup, lqd, method = 'natural')
    lqd_exp = function(t) exp(lqd_sp(t))
    # Get grid for density space
    dtemp = t0 + c(0, cumsum(sapply(2:length(lqdSup), function(i) integrate(lqd_exp, lqdSup[i - 1], lqdSup[i])$value))) - integrate(lqd_exp, lqdSup[1], lqdSup[c_ind])$value
    
  } else {
    # Get grid and function for density space
    dtemp = t0 + cumtrapzRcpp(lqdSup, exp(lqd)) - trapzRcpp(lqdSup[1:c_ind], exp(lqd[1:c_ind]))
  }
  
  # Remove duplicates
  indL = duplicated(dtemp[1:floor(M/2)], fromLast = TRUE)
  indR = duplicated(dtemp[(floor(M/2)+1):M])
  dtemp = dtemp[!c(indL, indR)]
  dens_temp = exp(-lqd[!c(indL, indR)]);
  
  # Interpolate to dSup and normalize
  dSup = seq(dtemp[1], dtemp[length(dtemp)], length.out = M)
  dens = approx(x = dtemp, y = dens_temp, xout = dSup, rule = c(2,2))[[2]]
  dens = dens/trapzRcpp(X = dSup,Y = dens)*(lqdSup[M] - lqdSup[1]); # Normalize, accounting for boundary cutoff
  
  return(list('dSup' = dSup, 'dens' = dens))
}

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
par(mfrow = c(1, 2))
matplot(u, dens, type = 'l', main = 'Original Density Estimates')
matplot(u, dens2, type = 'l', main = 'Density Estimates after Forward/Backward Transform')

par(mfrow = c(1, 1))
i = 10
plot(u, dens2[,i], type='l', col='red')
lines(u, dens[,i], type='l', col='black')




# SAVE DATAFRAME
i = 10
df_comparison <- data.frame(u = u, 
                     original_dens = dens[,i],
                     backwards_dens = dens2[,i])
install.packages("writexl")
library(writexl)
# Write the data frame to an Excel file
write_xlsx(x = df_comparison, path = "/Users/vizzotto/Documents/GitHub/densities4risk/data/processed/r_transformation.xlsx")



# BOVESPA
install.packages("readxl")
library(readxl)

bovespa_t <- read_excel("/Users/vizzotto/Documents/GitHub/densities4risk/data/interim/df_for_r.xlsx")
plot(bovespa_t$u, bovespa_t$density, type='l')
n = ncol(bovespa_t)
N = length(bovespa_t$u)

M = 5001
u <- bovespa_t$u
t0 = u[which.min(abs(u))]
dens_t <- bovespa_t$density
t = seq(0, 1, length.out = M)
lqdensity <- dens2lqd(dens = dens_t, dSup = u, lqdSup = t, t0 = t0, verbose = FALSE)
plot(t, lqdensity$lqd, type='l')
dens1 <- lqd2dens(lqd = lqdensity$lqd, lqdSup = t, t0 = t0, c = lqdensity$c, useSplines = FALSE, verbose = FALSE)
plot(dens1$dSup, dens1$dens, type='l', col='red')
dens2 <- approx(x = dens1$dSup, y = dens1$dens, xout = u, yleft = 0, yright = 0)$y
par(mfrow = c(1, 1))
plot(u, dens2, type='l', col='red')
lines(u, dens_t, type='l', col='black')
bovespa_t$backwards_density <- dens2
write_xlsx(x = bovespa_t, path = "/Users/vizzotto/Documents/GitHub/densities4risk/data/processed/r_transformation_bovespa_t.xlsx")
