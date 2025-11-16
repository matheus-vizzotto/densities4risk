
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


data(DJI_return)
library("ftsa")
library(fdapace) # -> trapzRcpp 
data <- DJI_return

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


### DADOS DA PESQUISA
library("readxl")
Y0 <- read_excel("C:/Users/user/Projetos/densities4risk/data/processed/kde.xlsx")
Y <- as.matrix(Y0[ , -1])



n = ncol(Y)
M = nrow(Y) # number of gridpoints for LQD functions - chosen large here so that 0 isn't too close to the boundary of all supports
u = seq(from = -0.01032280959030984, to = 0.01017648430576656, length = M)
lqd = matrix(0, nrow = M, ncol = n)
c = rep(0, 1, n)
t = seq(0, 1, length.out = M)
t0 = u[which.min(abs(u))] # closest value to 0
colnames(Y) <- 1:ncol(Y)

for(i in 1:n)
{
    tmp = dens2lqd(dens = Y[,i], dSup = u, lqdSup = t, t0 = t0, verbose = FALSE)
    lqd[,i] = tmp$lqd
    c[i] = tmp$c
}