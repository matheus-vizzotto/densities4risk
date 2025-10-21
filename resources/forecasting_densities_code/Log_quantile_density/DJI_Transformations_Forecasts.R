# Demonstrate modified LQD transformation for DJI cross-sectional return densities

setwd("/Volumes/MY PASSPORT/Dropbox/Todos/cs/code/DJI/Log_quantile_density")
source('dens2lqdNew.R')
source('lqd2densNew.R')

setwd("/Volumes/MY PASSPORT/Dropbox/Todos/cs/code/DJI/Horta-Ziegelman")
source("super_fun.R")

library(fdapace) # has trapezoidal integration used in the transformation functions, also FPCA
require(vars)
require(flexmix)
require(boot)
require(ftsa)

# Importing the data
setwd("/Volumes/MY PASSPORT/Dropbox/Todos/cs/data")
DJI_data = read.table("DJI_monthly.txt", header = TRUE, sep = "")
DJI_ret = as.matrix(DJI_data[,2:31])

Alex <- function(data, h_scale = 1, M = 3001,  m = 5001, lag_maximum = 4, no_boot = 1000, 
                 alpha_val = 0.1, p = 5, band_choice = c("Silverman", "DPI"), 
                 kernel = c("gaussian", "epanechnikov"), 
                 forecasting_method = c("uni", "BYZ"), var_prop = 0.85, fmethod)
{
    kernel = match.arg(kernel)
    forecasting_method = match.arg(forecasting_method)
    
    # Sample size
    N = nrow(data)
    
    if (!exists('h_scale')) h_scale = 1
    if(band_choice == "Silverman")
    {
        if(kernel == "gaussian")
        {
            h.hat_5m = sapply(1:N, function(t) 1.06*sd(data[t,])*(length(data[t,])^(-(1/5))))
        }
        if(kernel == "epanechnikov")
        {
            h.hat_5m = sapply(1:N, function(t) 2.34*sd(data[t,])*(length(data[t,])^(-(1/5))))  
        }
        h.hat_5m = h_scale * h.hat_5m
    }
    if(band_choice == "DPI")
    {
        if(kernel == "gaussian")
        {
            h.hat_5m = sapply(1:N, function(t) dpik(data[t,], kernel = "normal"))
        }
        if(kernel == "epanechnikov")
        {
            h.hat_5m = sapply(1:N, function(t) dpik(data[t,], kernel = "epanech"))
        }
        h.hat_5m = h_scale * h.hat_5m
    }
    
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
    dens = sapply(1:n, function(i) Y[,i]/trapzRcpp(X = u, Y = Y[,i]))

    # Try Forward transformation
    # number of gridpoints for LQD functions - chosen large here so that 0 isn't too close to the boundary of all supports
    lqd = matrix(0, nrow = M, ncol = n)
    c = vector("numeric", n)
    t = seq(0, 1, length.out = M)
    t0 = u[which.min(abs(u))] # closest value to 0
    for(i in 1:n)
    {
        tmp = dens2lqd(dens = dens[,i], dSup = u, lqdSup = t, t0 = t0, verbose = FALSE)
        lqd[,i] = tmp$lqd
        c[i] = tmp$c
    }
    
    ################################
    # Reconstruction densities
    # Now Try Backward Transform
    # Try cutting off boundary points with large LQD values (effectively setting the density to zero rather than trying to compute it numerically)
    ################################
    
    cut = res2 = list()
    for(i in 1:n)
    { 
        cut[[i]] = c(0, 0)
        cut[[i]][1] = sum(lqd[1:15,i] > 7)
        cut[[i]][2] = sum(lqd[(M-14):M, i] > 7)
        res2[[i]] = lqd2dens(lqd = lqd[,i], lqdSup = t, t0 = t0, c = c[i], cut = cut[[i]], useSplines = FALSE, verbose = FALSE)
    }
    dens2 = sapply(res2, function(r) approx(x = r$dSup, y = r$dens, xout = u, yleft = 0, yright = 0)$y)
  
    # Assess loss of mass incurred by boundary cutoff
    totalMass = apply(dens2, 2, function(d) trapzRcpp(X = u, Y = d))
    
    # Numerical comparison of densities 

    L2Diff = sapply(1:n, function(i) sqrt(trapzRcpp(X = u, Y = (dens[,i] - dens2[,i])^2))) # L^2 norm difference
    unifDiff = sapply(1:n, function(i){
                                interior = which(dens2[,i] > 0)
                                max(abs(dens[interior,i] - dens2[interior,i]))}) 
    # Uniform Metric excluding missing boundary values (due to boundary cutoff)
    
    ######################
    # Forecasting density
    ######################
    
    c_fore = as.numeric(inv.logit(forecast(auto.arima(logit(c)), h = 1)$mean))
    
    # cut-off two boundary points
    
    t_new = t[2:(M-1)]
    lqd_new = lqd[2:(M-1), ]
    
    if(forecasting_method == "uni")
    {
        ftsm_fitting = ftsm(y = fts(x = t_new, y = lqd_new), order = 10)
        ftsm_ncomp = head(which(cumsum((ftsm_fitting$varprop^2)/sum((ftsm_fitting$varprop^2))) >= var_prop),1)
        den_fore = forecast(object = ftsm(y = fts(x = t_new, y = lqd_new), order = ftsm_ncomp), h = 1, method = fmethod)$mean$y
    }
    if(forecasting_method == "BYZ")
    {
        dt = diff(t)[1]
        foo_out = super_fun(Y = lqd_new, lag_max = lag_maximum, B = no_boot, alpha = alpha_val, du = dt, 
                            p = p, m = (M-2), u = t_new, select_ncomp = "TRUE")
      
        # read outputs
      
        Ybar_est = foo_out$Ybar
        psihat_est = foo_out$psihat
        etahat_est = matrix(foo_out$etahat, ncol = n)
        selected_d0 = foo_out$d0
        
        score_object = t(etahat_est)
        colnames(score_object) = 1:selected_d0
        if(selected_d0 == 1)
        {
            etahat_pred_val = forecast(auto.arima(as.numeric(score_object)), h = 1)$mean
        }  
        else
        {
            VAR_mod = VARselect(score_object)
            etahat_pred = predict(VAR(y = score_object, p = min(VAR_mod$selection[3], 3), type = VAR_type), n.ahead = 1)
            etahat_pred_val = as.matrix(sapply(1:selected_d0, function(t) (etahat_pred$fcst[[t]])[1]))
        }
        den_fore = Ybar_est + psihat_est %*% etahat_pred_val
    }
    
    # add two useless boundary points back
    
    den_fore = as.matrix(c(8, den_fore, 8))
    cut = c(sum(den_fore[1:15,1] > 7), sum(den_fore[(M-14):M,1] > 7))
    
    res_fore = lqd2dens(lqd = den_fore, lqdSup = t, t0 = t0, c = c_fore, cut = cut, useSplines = FALSE, verbose = FALSE)
    dens_fore = approx(x = res_fore$dSup, y = res_fore$dens, xout = u, yleft = 0, yright = 0)$y
    return(list(L2Diff = L2Diff, unifDiff = unifDiff, density_reconstruct = dens2, density_original = dens,
                dens_fore = dens_fore, totalMass = range(totalMass), u = u))
}

dum = Alex(data = scale(US_female_pop, center = TRUE), forecasting_method = "uni", fmethod = "ets", kernel = "gaussian")


dum = Alex(data = DJI_ret, forecasting_method = "uni", fmethod = "ets", kernel = "gaussian")
dum = Alex(data = DJI_ret, forecasting_method = "uni", fmethod = "ets", kernel = "epanechnikov")

savepdf("Density_transform", width = 20, height = 10, toplines = 0.6)
par(mfrow = c(1,2))
plot(fts(dum$u, dum$density_original), xlab = "Grid points", ylab = "Density", main = "Original density")
plot(fts(dum$u, dum$density_reconstruct), xlab = "Grid points", ylab = "Density", main = "Reconstructed density")
dev.off()

# plot(fts(dum$u, dum$dens_fore), xlab = "Grid points", ylab = "Forecast density")


##########
# Holdout 
##########

ken_density = function(data, m = 5001, band_choice = c("Silverman", "DPI"), 
                       kernel = c("gaussian", "epanechnikov"))
{
    kernel = match.arg(kernel)
    # Sample size
    N = nrow(data)
    
    if (!exists('h_scale')) h_scale = 1
    if(band_choice == "Silverman")
    {
        if(kernel == "gaussian")
        {
            h.hat_5m = sapply(1:N, function(t) 1.06*sd(data[t,])*(length(data[t,])^(-(1/5))))
        }
        if(kernel == "epanechnikov")
        {
            h.hat_5m = sapply(1:N, function(t) 2.34*sd(data[t,])*(length(data[t,])^(-(1/5))))  
        }
        h.hat_5m = h_scale * h.hat_5m
    }
    if(band_choice == "DPI")
    {
        if(kernel == "gaussian")   
        {
            h.hat_5m = sapply(1:N, function(t) dpik(data[t,], kernel = "normal"))
        }
        if(kernel == "epanechnikov")
        {
            h.hat_5m = sapply(1:N, function(t) dpik(data[t,], kernel = "epanech"))
        }
        h.hat_5m = h_scale * h.hat_5m
    }
    
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
    return(list(Y = Y, u = u))
}

############
# Forecasts
############

test_density = ken_density(data = DJI_ret, band_choice = "Silverman", kernel = "gaussian")$Y[,111:165]
n_test = ncol(test_density)
test_density_epan = ken_density(data = DJI_ret, band_choice = "Silverman", kernel = "epanechnikov")$Y[,111:165]

test_density_DPI = ken_density(data = DJI_ret, band_choice = "DPI", kernel = "gaussian")$Y[,111:165]
test_density_epan_DPI = ken_density(data = DJI_ret, band_choice = "DPI", kernel = "epanechnikov")$Y[,111:165]


######################################
# Forecasts (univariate FTS approach)
######################################

# fmethod = "ets"

Alex_ets_DJI_ret_fore_den = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_ets_DJI_ret_fore_den[,ik] = Alex(data = DJI_ret[1:(109+ik),], band_choice = "Silverman", 
                                          forecastign_method = "uni", fmethod = "ets", 
                                          kernel = "gaussian")$dens_fore
    print(109+ik)
}  

Alex_ets_DJI_ret_fore_den_epan = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_ets_DJI_ret_fore_den_epan[,ik] = Alex(data = DJI_ret[1:(109+ik),], band_choice = "Silverman", 
                                               forecasting_method = "uni", fmethod = "ets", 
                                               kernel = "epanechnikov")$dens_fore
}

# fmethod = "arima"

Alex_arima_DJI_ret_fore_den = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_arima_DJI_ret_fore_den[,ik] = Alex(data = DJI_ret[1:(109+ik),], band_choice = "Silverman", 
                                            forecasting_method = "uni", fmethod = "arima", 
                                            kernel = "gaussian")$dens_fore
    print(109+ik)
}  

Alex_arima_DJI_ret_fore_den_epan = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_arima_DJI_ret_fore_den_epan[,ik] = Alex(data = DJI_ret[1:(109+ik),], band_choice = "Silverman", 
                                                 forecasting_method = "uni", fmethod = "arima", 
                                                 kernel = "epanechnikov")$dens_fore
}


###########################
# Forecasts (BYZ approach)
###########################

Alex_BYZ_DJI_ret_fore_den = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_BYZ_DJI_ret_fore_den[,ik] = Alex(data = DJI_ret[1:(109+ik),], band_choice = "Silverman", 
                                          forecasting_method = "BYZ", kernel = "gaussian")$dens_fore
    print(109+ik)
}

Alex_BYZ_DJI_ret_fore_den_epan = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_BYZ_DJI_ret_fore_den_epan[,ik] = Alex(data = DJI_ret[1:(109+ik),], band_choice = "Silverman", 
                                               forecasting_method = "BYZ", kernel = "epanechnikov")$dens_fore
    print(109+ik)
}

##########
# Holdout 
##########

ken_density = function(data, m = 5001, band_choice = c("Silverman", "DPI"), 
                       kernel = c("gaussian", "epanechnikov"))
{
    kernel = match.arg(kernel)
    # Sample size
    N = nrow(data)

    if (!exists('h_scale')) h_scale = 1
    if(band_choice == "Silverman")
    {
        if(kernel == "gaussian")
        {
            h.hat_5m = sapply(1:N, function(t) 1.06*sd(data[t,])*(length(data[t,])^(-(1/5))))
        }
        if(kernel == "epanechnikov")
        {
            h.hat_5m = sapply(1:N, function(t) 2.34*sd(ret[[t]])*(length(ret[[t]])^(-(1/5))))  
        }
        h.hat_5m = h_scale * h.hat_5m
    }
    if(band_choice == "DPI")
    {
        if(kernel == "gaussian")
        {
            h.hat_5m = sapply(1:N, function(t) dpik(data[t,], kernel = "normal"))
        }
        if(kernel == "epanechnikov")
        {
            h.hat_5m = sapply(1:N, function(t) dpik(data[t,], kernel = "epanech"))
        }
        h.hat_5m = h_scale * h.hat_5m
    }
    
    
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
    return(list(Y = Y, u = u))
}

############
# Forecasts
############

test_density = ken_density(data = DJI_ret, band_choice = "Silverman", kernel = "gaussian")$Y[,111:165]
n_test = ncol(test_density)
test_density_epan = ken_density(data = DJI_ret, band_choice = "Silverman", kernel = "epanechnikov")$Y[,111:165]

test_density_DPI = ken_density(data = DJI_ret, band_choice = "DPI", kernel = "gaussian")$Y[,111:165]
test_density_epan_DPI = ken_density(data = DJI_ret, band_choice = "DPI", kernel = "epanechnikov")$Y[,111:165]


##############################
# Kullback-Leibler divergence
##############################

# kernel = "gaussian"

KLdiv_Alex_ets = KLdiv_Alex_arima = KLdiv_Alex_BYZ = matrix(NA, n_test, 2)
for(ik in 1:55)
{
    Alex_ets_compar = cbind(test_density[,ik], Alex_ets_DJI_ret_fore_den[,ik])
    Alex_arima_compar = cbind(test_density[,ik], Alex_arima_DJI_ret_fore_den[,ik])
    Alex_BYZ_compar = cbind(test_density[,ik], Alex_BYZ_DJI_ret_fore_den[,ik])
    colnames(Alex_ets_compar) = colnames(Alex_arima_compar) = colnames(Alex_BYZ_compar) = c("True", "Estimate")
    
    KLdiv_Alex_ets[ik,] = as.numeric(KLdiv(Alex_ets_compar, eps = 1e-16))[2:3]
    KLdiv_Alex_arima[ik,] = as.numeric(KLdiv(Alex_arima_compar, eps = 1e-16))[2:3]
    KLdiv_Alex_BYZ[ik,] = as.numeric(KLdiv(Alex_BYZ_compar, eps = 1e-16))[2:3]
    print(ik)
    rm(Alex_ets_compar); rm(Alex_arima_compar); rm(Alex_BYZ_compar)
}

round(colMeans(KLdiv_Alex_ets), 4) # 0.3389 0.7206
round(colMeans(KLdiv_Alex_arima), 4) # 0.3392 0.7029
round(colMeans(KLdiv_Alex_BYZ), 4) # 0.3391 0.7883

round(sum(colMeans(KLdiv_Alex_ets)), 4) # 1.0594
KLdiv_Alex_arima_summary = round(sum(colMeans(KLdiv_Alex_arima)), 4) # 1.0421
round(sum(colMeans(KLdiv_Alex_BYZ)), 4) # 1.1273

# kernel = "epan"

KLdiv_Alex_ets_epan = KLdiv_Alex_arima_epan = KLdiv_Alex_BYZ_epan = matrix(NA, n_test, 2)
for(ik in 1:55)
{
    Alex_ets_compar = cbind(test_density_epan[,ik], Alex_ets_DJI_ret_fore_den_epan[,ik])
    Alex_arima_compar = cbind(test_density_epan[,ik], Alex_arima_DJI_ret_fore_den_epan[,ik])
    Alex_BYZ_compar = cbind(test_density_epan[,ik], Alex_BYZ_DJI_ret_fore_den_epan[,ik])
    colnames(Alex_ets_compar) = colnames(Alex_arima_compar) = colnames(Alex_BYZ_compar) = c("True", "Estimate")
    
    KLdiv_Alex_ets_epan[ik,] = as.numeric(KLdiv(Alex_ets_compar, eps = 1e-16))[2:3]
    KLdiv_Alex_arima_epan[ik,] = as.numeric(KLdiv(Alex_arima_compar, eps = 1e-16))[2:3]
    KLdiv_Alex_BYZ_epan[ik,] = as.numeric(KLdiv(Alex_BYZ_compar, eps = 1e-16))[2:3]
    print(ik)
    rm(Alex_ets_compar); rm(Alex_arima_compar); rm(Alex_BYZ_compar)
}

round(sum(colMeans(KLdiv_Alex_ets_epan)), 4) # 1.1776
KLdiv_Alex_arima_summary_epan = round(sum(colMeans(KLdiv_Alex_arima_epan)), 4) # 1.1833
round(sum(colMeans(KLdiv_Alex_BYZ_epan)), 4) # 1.2211

############################
# Jensen-Shannon divergence
############################

# kernel = "gaussian"

JSdiv_Alex_ets = JSdiv_Alex_arima = JSdiv_Alex_BYZ = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Alex_ets_compar = cbind(test_density[,ik], Alex_ets_DJI_ret_fore_den[,ik])
    M = rowMeans(Alex_ets_compar)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(Alex_ets_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_ets[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_ets_compar)
    
    Alex_arima_compar = cbind(test_density[,ik], Alex_arima_DJI_ret_fore_den[,ik])
    M = rowMeans(Alex_arima_compar)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(Alex_arima_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_arima[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_arima_compar)
    
    Alex_BYZ_compar = cbind(test_density[,ik], Alex_BYZ_DJI_ret_fore_den[,ik])
    M = rowMeans(Alex_BYZ_compar)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(Alex_BYZ_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_BYZ[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_BYZ_compar)
    print(ik)
}

round(sum(JSdiv_Alex_ets), 4) # 3.0353
JSdiv_Alex_arima_summary = round(sum(JSdiv_Alex_arima), 4) # 3.0129
round(sum(JSdiv_Alex_BYZ), 4) # 3.0197

# kernel = "epan"

JSdiv_Alex_ets_epan = JSdiv_Alex_arima_epan = JSdiv_Alex_BYZ_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Alex_ets_compar = cbind(test_density_epan[,ik], Alex_ets_DJI_ret_fore_den_epan[,ik])
    M = rowMeans(Alex_ets_compar)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(Alex_ets_DJI_ret_fore_den_epan[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_ets_epan[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_ets_compar)
    
    Alex_arima_compar = cbind(test_density_epan[,ik], Alex_arima_DJI_ret_fore_den_epan[,ik])
    M = rowMeans(Alex_arima_compar)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(Alex_arima_DJI_ret_fore_den_epan[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_arima_epan[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_arima_compar)
    
    Alex_BYZ_compar = cbind(test_density_epan[,ik], Alex_BYZ_DJI_ret_fore_den_epan[,ik])
    M = rowMeans(Alex_BYZ_compar)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(Alex_BYZ_DJI_ret_fore_den_epan[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_BYZ_epan[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_BYZ_compar)
    print(ik)
}

round(sum(JSdiv_Alex_ets_epan), 4) # 2.0183
JSdiv_Alex_arima_summary_epan = round(sum(JSdiv_Alex_arima_epan), 4) # 2.0165
round(sum(JSdiv_Alex_BYZ_epan), 4) # 2.0211

##################################
# Jensen-Shannon divergence (geo)
##################################

# kernel = "gaussian"

JSdiv_geo_Alex_ets = JSdiv_geo_Alex_arima = JSdiv_geo_Alex_BYZ = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Alex_ets_compar = cbind(test_density[,ik], Alex_ets_DJI_ret_fore_den[,ik])
    M = apply(Alex_ets_compar, 1, geometric.mean)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(Alex_ets_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_ets[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_ets_compar)
    
    Alex_arima_compar = cbind(test_density[,ik], Alex_arima_DJI_ret_fore_den[,ik])
    M = apply(Alex_arima_compar, 1, geometric.mean)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(Alex_arima_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_arima[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_arima_compar)
    
    Alex_BYZ_compar = cbind(test_density[,ik], Alex_BYZ_DJI_ret_fore_den[,ik])
    M = apply(Alex_BYZ_compar, 1, geometric.mean)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(Alex_BYZ_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_BYZ[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_BYZ_compar)
    print(ik)
}

JSdiv_geo_Alex_arima_summary = round(sum(JSdiv_geo_Alex_arima), 4) # 6.9443
JSdiv_geo_Alex_ets_summary = round(sum(JSdiv_geo_Alex_ets), 4)     # 7.0284
JSdiv_geo_Alex_BYZ_summary = round(sum(JSdiv_geo_Alex_BYZ), 4)     # 7.4372

# kernel = "epan"

JSdiv_geo_Alex_ets_epan = JSdiv_geo_Alex_arima_epan = JSdiv_geo_Alex_BYZ_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Alex_ets_compar = cbind(test_density_epan[,ik], Alex_ets_DJI_ret_fore_den_epan[,ik])
    M = apply(Alex_ets_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(Alex_ets_DJI_ret_fore_den_epan[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_ets_epan[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3] 		
    rm(M); rm(P_M); rm(E_M); rm(Alex_ets_compar)
    
    Alex_arima_compar = cbind(test_density_epan[,ik], Alex_arima_DJI_ret_fore_den_epan[,ik])
    M = apply(Alex_arima_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(Alex_arima_DJI_ret_fore_den_epan[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_arima_epan[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_arima_compar)
    
    Alex_BYZ_compar = cbind(test_density_epan[,ik], Alex_BYZ_DJI_ret_fore_den_epan[,ik])
    M = apply(Alex_BYZ_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(Alex_BYZ_DJI_ret_fore_den_epan[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_BYZ_epan[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_BYZ_compar)
    print(ik)
}

JSdiv_geo_Alex_arima_epan_summary = round(sum(JSdiv_geo_Alex_arima_epan), 4) # 7.0648
JSdiv_geo_Alex_ets_epan_summary   = round(sum(JSdiv_geo_Alex_ets_epan), 4) # 7.0322
JSdiv_geo_Alex_BYZ_epan_summary   = round(sum(JSdiv_geo_Alex_BYZ_epan), 4) # 7.2825


############
# L1 - norm
############

# kernel = "gaussian"

L1_norm_Alex_ets = L1_norm_Alex_arima = L1_norm_Alex_BYZ = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_Alex_ets[ik] = sum(abs(Alex_ets_DJI_ret_fore_den[,ik] - test_density[,ik]))
    L1_norm_Alex_arima[ik] = sum(abs(Alex_arima_DJI_ret_fore_den[,ik] - test_density[,ik]))
    L1_norm_Alex_BYZ[ik] = sum(abs(Alex_BYZ_DJI_ret_fore_den[,ik] - test_density[,ik]))
}  

round(mean(L1_norm_Alex_ets), 4) # 955.2805
L1_norm_Alex_arima_summary = round(mean(L1_norm_Alex_arima), 4) # 948.7655
round(mean(L1_norm_Alex_BYZ), 4) # 947.7571

# kernel = "epan"

L1_norm_Alex_ets_epan = L1_norm_Alex_arima_epan = L1_norm_Alex_BYZ_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_Alex_ets_epan[ik] = sum(abs(Alex_ets_DJI_ret_fore_den_epan[,ik] - test_density_epan[,ik]))
    L1_norm_Alex_arima_epan[ik] = sum(abs(Alex_arima_DJI_ret_fore_den_epan[,ik] - test_density_epan[,ik]))
    L1_norm_Alex_BYZ_epan[ik] = sum(abs(Alex_BYZ_DJI_ret_fore_den_epan[,ik] - test_density_epan[,ik]))
}

round(mean(L1_norm_Alex_ets_epan), 4) # 721.5399
L1_norm_Alex_arima_summary_epan = round(mean(L1_norm_Alex_arima_epan), 4) # 720.0048
round(mean(L1_norm_Alex_BYZ_epan), 4) # 720.0772

############
# L2 - norm
############

# kernel = "gaussian"

L2_norm_Alex_ets = L2_norm_Alex_arima = L2_norm_Alex_BYZ = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_Alex_ets[ik] = sqrt(sum((test_density[,ik] - Alex_ets_DJI_ret_fore_den[,ik])^2))
    L2_norm_Alex_arima[ik] = sqrt(sum((test_density[,ik] - Alex_arima_DJI_ret_fore_den[,ik])^2))
    L2_norm_Alex_BYZ[ik] = sqrt(sum((test_density[,ik] - Alex_BYZ_DJI_ret_fore_den[,ik])^2))
}

round(mean(L2_norm_Alex_ets), 4) # 46.2625
L2_norm_Alex_arima_summary = round(mean(L2_norm_Alex_arima), 4) # 45.8322
round(mean(L2_norm_Alex_BYZ), 4) # 45.7541

# kernel = "epan"

L2_norm_Alex_ets_epan = L2_norm_Alex_arima_epan = L2_norm_Alex_BYZ_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_Alex_ets_epan[ik] = sqrt(sum((test_density_epan[,ik] - Alex_ets_DJI_ret_fore_den_epan[,ik])^2))
    L2_norm_Alex_arima_epan[ik] = sqrt(sum((test_density_epan[,ik] - Alex_arima_DJI_ret_fore_den_epan[,ik])^2))
    L2_norm_Alex_BYZ_epan[ik] = sqrt(sum((test_density_epan[,ik] - Alex_BYZ_DJI_ret_fore_den_epan[,ik])^2))
}

round(mean(L2_norm_Alex_ets_epan), 4) # 29.1429
L2_norm_Alex_arima_summary_epan = round(mean(L2_norm_Alex_arima_epan), 4) # 29.0229
round(mean(L2_norm_Alex_BYZ_epan), 4) # 29.0252

##############
# Linf - norm
##############

# kernel = "gaussian"

Linf_norm_Alex_ets = Linf_norm_Alex_arima = Linf_norm_Alex_BYZ = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_Alex_ets[ik] = max(abs(Alex_ets_DJI_ret_fore_den[,ik] - test_density[,ik]))
    Linf_norm_Alex_arima[ik] = max(abs(Alex_arima_DJI_ret_fore_den[,ik] - test_density[,ik]))
    Linf_norm_Alex_BYZ[ik] = max(abs(Alex_BYZ_DJI_ret_fore_den[,ik] - test_density[,ik]))
}  

round(mean(Linf_norm_Alex_ets), 4) # 3.6844
Linf_norm_Alex_arima_summary = round(mean(Linf_norm_Alex_arima), 4) # 3.6546  
round(mean(Linf_norm_Alex_BYZ), 4) # 3.6469

# kernel = "epan"

Linf_norm_Alex_ets_epan = Linf_norm_Alex_arima_epan = Linf_norm_Alex_BYZ_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_Alex_ets_epan[ik] = max(abs(Alex_ets_DJI_ret_fore_den_epan[,ik] - test_density_epan[,ik]))
    Linf_norm_Alex_arima_epan[ik] = max(abs(Alex_arima_DJI_ret_fore_den_epan[,ik] - test_density_epan[,ik]))
    Linf_norm_Alex_BYZ_epan[ik] = max(abs(Alex_BYZ_DJI_ret_fore_den_epan[,ik] - test_density_epan[,ik]))
}

round(mean(Linf_norm_Alex_ets_epan), 4)  # 1.7815
Linf_norm_Alex_arima_summary_epan = round(mean(Linf_norm_Alex_arima_epan), 4)  # 1.772
round(mean(Linf_norm_Alex_BYZ_epan), 4)  # 1.7731

##########
# summary
##########

DJI_LQD_summary = c(KLdiv_Alex_arima_summary,
                    JSdiv_Alex_arima_summary,
                    JSdiv_geo_Alex_arima_summary,
                    L1_norm_Alex_arima_summary,
                    L2_norm_Alex_arima_summary,
                    Linf_norm_Alex_arima_summary)

DJI_LQD_summary_epan = c(KLdiv_Alex_arima_summary_epan,
                         JSdiv_Alex_arima_summary_epan,
                         JSdiv_geo_Alex_arima_epan_summary,
                         L1_norm_Alex_arima_summary_epan,
                         L2_norm_Alex_arima_summary_epan,
                         Linf_norm_Alex_arima_summary_epan)

