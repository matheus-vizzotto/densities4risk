# load R package

require(vars)

# Setting up the functions

# Defining the inner product on L^2([a,b])
inner.product = function(f, g, du) {drop(crossprod(f, g))*du}

# Defining the L2-norm
L2norm = function(f, du) {sqrt(inner.product(f, f, du))}

# Importing the data
setwd("/Volumes/MY PASSPORT/Dropbox/Todos/cs/data")
DJI_data = read.table("DJI_monthly.txt", header = TRUE, sep = "")
DJI_ret = as.matrix(DJI_data[,2:31])

DJI_date = DJI_data[,1]

setwd("/Volumes/MY PASSPORT/Dropbox/Todos/cs/code/DJI/Horta-Ziegelman")
source("super_fun.R")

# data: DJI_ret
# h_scale: scaling parameter in the kernel density estimator
# p: the number of backward parameters
# m: the number of grid points
# D: the number of retained components

Horta_Ziegelman <- function(data, h_scale = 1, p = 5, m = 5001, lag_maximum = 6, no_boot = 1000,
                            alpha_val = 0.10, kernel = c("gaussian", "epanechnikov"), 
                            band_choice = c("Silverman", "DPI"),
                            VAR_type = "both", ncomp_select = "TRUE", D_val = 10)
{
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
    
    # Initialization parameters
    n = N # Number of daily observations
    
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
    
    # main function
    
    foo_out = super_fun(Y = Y, lag_max = lag_maximum, B = no_boot, alpha = alpha_val, du = du, 
                        p = p, m = m, u = u, select_ncomp = ncomp_select, dimension = D_val)
    
    # read outputs
    
    Ybar_est = foo_out$Ybar
    psihat_est = foo_out$psihat
    etahat_est = matrix(foo_out$etahat, ncol = N)
    selected_d0 = foo_out$d0
    thetahat_val = foo_out$thetahat
    if(ncomp_select == "TRUE")
    {
        selected_d0_pvalues = foo_out$bs.pvalues
    }
    else
    {
        selected_d0_pvalues = 10^5 
    }
    
    # VAR forecasting
    
    require(vars)
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
    Yhat.fix_den = den_fore = Ybar_est + psihat_est %*% etahat_pred_val
    
    # adjustment
    
    Yhat.fix_den[Yhat.fix_den < 0] = 0
    Yhat.fix_den = Yhat.fix_den/(sum(Yhat.fix_den)*du)
    return(list(Yhat.fix_den = Yhat.fix_den, u = u, du = du, Ybar_est = Ybar_est,
                psihat_est = psihat_est, etahat_est = etahat_est, etahat_pred_val = etahat_pred_val, 
                selected_d0 = selected_d0, selected_d0_pvalues = selected_d0_pvalues, thetahat_val = thetahat_val))
}

HZ_den_fore = Horta_Ziegelman(data = DJI_ret, ncomp_select = "FALSE")

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
    return(list(Y = Y, u = u, du = du))
}

############
# Test set
# Forecasts
############

dum = ken_density(data = DJI_ret, band_choice = "Silverman", kernel = "gaussian")
DJI_test_density_grid = dum$u
test_density = dum$Y[,111:165]

savepdf("DJI_fts_plot", width = 12, height = 10, toplines = 1)
plot(fts(DJI_test_density_grid[3000:4501], 
      test_density[3000:4501,44:55] * dum$du),
      xlab = "Grid point", ylab = "Density", main = "Dow-Jones Index (January to December, 2017)")
dev.off()

n_test = ncol(test_density)

test_density_epan = ken_density(data = DJI_ret, band_choice = "Silverman", kernel = "epanechnikov")$Y[,111:165]

####################
# Density forecasts
####################

# Gaussian kernel

DJI_Horta_Ziegelman_ncomp_gaussian = c(3, 3, 3, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                       5, 1, 5, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 1, 1, 3, 
                                       3, 3, 3, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 3, 3, 4, 4, 4, 3)

DJI_ret_fore_den = matrix(NA, 5001, n_test)
DJI_ret_fore_u = DJI_ret_fore_d0 = vector("numeric", n_test)
DJI_ret_fore_d0_pvalues = matrix(NA, 6, n_test)
for(ik in 1:55)
{
    dum = Horta_Ziegelman(data = DJI_ret[1:(109+ik),], band_choice = "Silverman", kernel = "gaussian", ncomp_select = "FALSE", D_val = DJI_Horta_Ziegelman_ncomp_gaussian[ik])
    DJI_ret_fore_den[,ik] = dum$Yhat.fix_den
    DJI_ret_fore_u[ik] = dum$du
    DJI_ret_fore_d0[ik] = dum$selected_d0
    DJI_ret_fore_d0_pvalues[,ik] = dum$selected_d0_pvalues
    print(109+ik); rm(dum); rm(ik)
}  

# Epanechnikov kernel

DJI_Horta_Ziegelman_ncomp_epan = c(3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                   1, 3, 3, 3, 3, 3, 3, 3, 1, 3, 1, 3, 3, 3, 6, 6, 6, 6, 6, 6, 
                                   1, 6, 6, 6, 6, 6, 6, 6, 1, 1, 6, 1, 1, 1, 1)

DJI_ret_fore_den_epan = matrix(NA, 5001, n_test)
DJI_ret_fore_u_epan = DJI_ret_fore_d0_epan = vector("numeric", n_test)
DJI_ret_fore_d0_pvalues_epan = matrix(NA, 6, n_test)
for(ik in 1:n_test)
{
    dum = Horta_Ziegelman(data = DJI_ret[1:(109+ik),], band_choice = "Silverman", kernel = "epanechnikov", ncomp_select = "FALSE", D_val = DJI_Horta_Ziegelman_ncomp_epan[ik])
    DJI_ret_fore_den_epan[,ik] = dum$Yhat.fix_den
    DJI_ret_fore_u_epan[ik] = dum$du
    DJI_ret_fore_d0_epan[ik] = dum$selected_d0
    DJI_ret_fore_d0_pvalues_epan[,ik] = dum$selected_d0_pvalues
    print(109+ik); rm(dum); rm(ik)
}

#####################
# Density evaluation
#####################

## Kullback-Leibler divergence

# kernel = "gaussian"

KLdiv_Horta_Ziegelman = matrix(NA, n_test, 2)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density[,ik], DJI_ret_fore_den[,ik])
    colnames(Horta_Ziegelman_compar) = c("True", "Estimate")
    KLdiv_Horta_Ziegelman[ik,] = as.numeric(KLdiv(Horta_Ziegelman_compar, eps=1e-16))[2:3]
    print(ik)
    rm(Horta_Ziegelman_compar)
}

round(colMeans(KLdiv_Horta_Ziegelman), 4) 
KLdiv_Horta_Ziegelman_summary = round(sum(colMeans(KLdiv_Horta_Ziegelman)), 4) # 1.3070

# kernel = "epan"

KLdiv_Horta_Ziegelman_epan = matrix(NA, n_test, 2)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density_epan[,ik], DJI_ret_fore_den_epan[,ik])
    colnames(Horta_Ziegelman_compar) = c("True", "Estimate")
    KLdiv_Horta_Ziegelman_epan[ik,] = as.numeric(KLdiv(Horta_Ziegelman_compar, eps=1e-16))[2:3]
    print(ik)
    rm(Horta_Ziegelman_compar)
}

round(colMeans(KLdiv_Horta_Ziegelman_epan), 4) # 1.2252 0.2473
KLdiv_Horta_Ziegelman_summary_epan = round(sum(colMeans(KLdiv_Horta_Ziegelman_epan)), 4) # 1.4725


## Jensen-Shannon divergence

# kernel = "gaussian"

JSdiv_Horta_Ziegelman = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density[,ik], DJI_ret_fore_den[,ik])
    M = rowMeans(Horta_Ziegelman_compar)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Horta_Ziegelman[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
}

JSdiv_Horta_Ziegelman_summary = round(sum(JSdiv_Horta_Ziegelman), 4) #  3.5986

# kernel = "epan"

JSdiv_Horta_Ziegelman_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density_epan[,ik], DJI_ret_fore_den_epan[,ik])
    M = rowMeans(Horta_Ziegelman_compar)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(DJI_ret_fore_den_epan[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Horta_Ziegelman_epan[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
}

JSdiv_Horta_Ziegelman_summary_epan = round(sum(JSdiv_Horta_Ziegelman_epan), 4) # 2.2652

## Jensen-Shannon divergence (geometric mean)

# kernel = "gaussian"

JSdiv_geo_Horta_Ziegelman = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density[,ik], DJI_ret_fore_den[,ik])
    M = apply(Horta_Ziegelman_compar, 1, geometric.mean)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Horta_Ziegelman[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    print(ik); rm(Horta_Ziegelman_compar)
}

JSdiv_geo_Horta_Ziegelman_summary = round(sum(JSdiv_geo_Horta_Ziegelman), 4) # 9.4038

# kernel = "epan"

JSdiv_geo_Horta_Ziegelman_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density_epan[,ik], DJI_ret_fore_den_epan[,ik])
    M = apply(Horta_Ziegelman_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(DJI_ret_fore_den_epan[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Horta_Ziegelman_epan[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    print(ik); rm(Horta_Ziegelman_compar)
}

JSdiv_geo_Horta_Ziegelman_summary_epan = round(sum(JSdiv_geo_Horta_Ziegelman_epan), 4) # 8.3987

## L1-norm

# kernel = "gaussian"

L1_norm_Horta_Ziegelman = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_Horta_Ziegelman[ik] = sum(abs(test_density[,ik] - DJI_ret_fore_den[,ik]))
}

L1_norm_Horta_Ziegelman_summary = round(mean(L1_norm_Horta_Ziegelman), 4) # 1039.357

# kernel = "epan"

L1_norm_Horta_Ziegelman_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_Horta_Ziegelman_epan[ik] = sum(abs(test_density_epan[,ik] - DJI_ret_fore_den_epan[,ik]))
}

L1_norm_Horta_Ziegelman_summary_epan = round(mean(L1_norm_Horta_Ziegelman_epan), 4) # 756.9432  

## L2-norm

# kernel = "gaussian"

L2_norm_Horta_Ziegelman = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_Horta_Ziegelman[ik] = sqrt(sum((test_density[,ik] - DJI_ret_fore_den[,ik])^2))
}

L2_norm_Horta_Ziegelman_summary = round(mean(L2_norm_Horta_Ziegelman), 4) # 46.9608

# kernel = "epan"

L2_norm_Horta_Ziegelman_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_Horta_Ziegelman_epan[ik] = sqrt(sum((test_density_epan[,ik] - DJI_ret_fore_den_epan[,ik])^2))
}

L2_norm_Horta_Ziegelman_summary_epan = round(mean(L2_norm_Horta_Ziegelman_epan), 4) # 29.2258

## Linf-norm

# kernel = "gaussian"

Linf_norm_Horta_Ziegelman = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_Horta_Ziegelman[ik] = max(abs(test_density[,ik] - DJI_ret_fore_den[,ik]))
}

Linf_norm_Horta_Ziegelman_summary = round(mean(Linf_norm_Horta_Ziegelman), 4) # 3.8713

# kernel = "epan"

Linf_norm_Horta_Ziegelman_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_Horta_Ziegelman_epan[ik] = max(abs(test_density_epan[,ik] - DJI_ret_fore_den_epan[,ik]))
}

Linf_norm_Horta_Ziegelman_summary_epan = round(mean(Linf_norm_Horta_Ziegelman_epan), 4) # 1.7911

###########
## summary
###########

DJI_Horta_Ziegelman_summary = c(KLdiv_Horta_Ziegelman_summary,
                                JSdiv_Horta_Ziegelman_summary,
                                JSdiv_geo_Horta_Ziegelman_summary,
                                L1_norm_Horta_Ziegelman_summary,
                                L2_norm_Horta_Ziegelman_summary,
                                Linf_norm_Horta_Ziegelman_summary)

DJI_Horta_Ziegelman_summary_epan = c(KLdiv_Horta_Ziegelman_summary_epan,
                                     JSdiv_Horta_Ziegelman_summary_epan,
                                     JSdiv_geo_Horta_Ziegelman_summary_epan,
                                     L1_norm_Horta_Ziegelman_summary_epan,
                                     L2_norm_Horta_Ziegelman_summary_epan,
                                     Linf_norm_Horta_Ziegelman_summary_epan)

