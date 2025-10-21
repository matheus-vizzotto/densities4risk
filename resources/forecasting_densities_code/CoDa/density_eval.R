########################
# Forecasting densities
########################

setwd("/Volumes/MY PASSPORT/Dropbox/Todos/cs/code/CoDa")
source("CoDa.R")
require(flexmix)
require(ftsa)
require(psych)

# Defining the inner product on L^2([a,b])
inner.product = function(f, g, du)
{
    return(drop(crossprod(f, g)) * du)
}

# Defining the L2-norm
L2norm = function(f, du) 
{
    return(sqrt(inner.product(f, f, du)))
}

# Importing the data
setwd("/Volumes/MY PASSPORT/Dropbox/Todos/cs/data")
DJI_data = read.table(file = "DJI_monthly.txt", header = TRUE, sep = "")
DJI_ret = as.matrix(DJI_data[,2:31])


DJI_forecast_den <- function(data, normalization, ik, h_scale = 1, m = 5001, 
                             band_choice = c("Silverman", "DPI"),
                             kernel = c("gaussian", "epanechnikov"), 
                             var_prop = 0.99)   
{
    band_choice = match.arg(band_choice)
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
    
    ###########################
    # Dealing with zero values
    ###########################
    
    DJI_return_density_train_trans <- Y[,1:(109+ik)]
    DJI_return_density_train_transformation = DJI_return_density_train_trans * (10^6)
    n_1 = ncol(DJI_return_density_train_transformation)
    epsilon = sapply(1:n_1, function(X) max(DJI_return_density_train_transformation[,X] - round(DJI_return_density_train_transformation[,X], 2)))
    
    DJI_CoDa_mat = matrix(NA, m, n_1)
    for(ik in 1:n_1)
    { 
        index = which(round(DJI_return_density_train_transformation[,ik], 2) == 0)
        DJI_CoDa_mat[,ik] = replace(DJI_return_density_train_transformation[,ik], index, epsilon[ik])
        DJI_CoDa_mat[-index,ik] = DJI_return_density_train_transformation[-index,ik] * (1 - (length(index) * epsilon[ik])/(10^6))
    }    
    
    # CoDa
    
    c = colSums(Y)[1]
    dum = CoDa_recon(dat = t(DJI_CoDa_mat), normalize = normalization, 
                     fore_method = "ETS", fh = 1, varprop = var_prop, constant = c)
    return(dum$d_x_t_star_fore)
}

################################
# obtain 55 estimated densities
################################

### normalization = TRUE

## band_choice = "Silverman"

# kernel = "gaussian"

require(doMC)
registerDoMC(8)
den_val = foreach(iw = 1:55) %dopar% DJI_forecast_den(data = DJI_ret, normalization = TRUE, 
                                                      band_choice = "Silverman", 
                                                      ik = iw, kernel = "gaussian")

den_fore = matrix(NA, 5001, 55)
for(ik in 1:55)
{
    den_fore[,ik] = den_val[[ik]]
}
rm(den_val)

# kernel = "epanechnikov"

registerDoMC(8)
den_val_epan = foreach(iw = 1:55) %dopar% DJI_forecast_den(data = DJI_ret, normalization = TRUE, 
                                                           band_choice = "Silverman", 
                                                           ik = iw, kernel = "epanechnikov") 

den_fore_epan = matrix(NA, 5001, 55)
for(ik in 1:55)
{
    den_fore_epan[,ik] = den_val_epan[[ik]]
}
rm(den_val_epan)

## band_choice = "DPI"

# kernel = "gaussian"

registerDoMC(8)
den_val_DPI = foreach(iw = 1:55) %dopar% DJI_forecast_den(data = DJI_ret, normalization = TRUE, 
                                                          band_choice = "DPI", 
                                                          ik = iw, kernel = "gaussian")

den_fore_DPI = matrix(NA, 5001, 55)
for(ik in 1:55)
{
    den_fore_DPI[,ik] = den_val_DPI[[ik]]
}
rm(den_val_DPI)

# kernel = "epanechnikov"

registerDoMC(8)
den_val_epan_DPI = foreach(iw = 1:55) %dopar% DJI_forecast_den(data = DJI_ret, normalization = TRUE,
                                                               band_choice = "DPI",
                                                               ik = iw, kernel = "epanechnikov")

den_fore_epan_DPI = matrix(NA, 5001, 55)
for(ik in 1:55)
{
    den_fore_epan_DPI[,ik] = den_val_epan_DPI[[ik]]
}
rm(den_val_epan_DPI)

#############################################################
# Keep the last 55 logspline densities as the testing sample
#############################################################

test_density = ken_density(data = DJI_ret, band_choice = "Silverman", kernel = "gaussian")$Y[,111:165]
test_density_epan = ken_density(data = DJI_ret, band_choice = "Silverman", kernel = "epanechnikov")$Y[,111:165]
n_test = ncol(test_density_epan)

test_density_DPI = ken_density(data = DJI_ret, band_choice = "DPI", kernel = "gaussian")$Y[,111:165]
test_density_epan_DPI = ken_density(data = DJI_ret, band_choice = "DPI", kernel = "epanechnikov")$Y[,111:165]

##############################
# Kullback-Leibler divergence
##############################

## band_choice = "Silverman"

# kernel = "gaussian"

KLdiv_CoDa = matrix(NA, 55, 2)
for(ik in 1:55)
{
    dat = cbind(true = test_density[,ik], forecast = den_fore[,ik])
    colnames(dat) = c("True", "Estimate")
    KLdiv_CoDa[ik,] = as.numeric(KLdiv(dat, eps=1e-16))[2:3]
    print(ik)
}

KLdiv_CoDa_summary = round(sum(colMeans(KLdiv_CoDa)), 4) # 0.6658

# kernel = "epan"

KLdiv_CoDa_epan = matrix(NA, 55, 2)
for(ik in 1:55)
{
    dat = cbind(true = test_density_epan[,ik], forecast = den_fore_epan[,ik])
    colnames(dat) = c("True", "Estimate")
    KLdiv_CoDa_epan[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
    print(ik)
    rm(dat)
}

KLdiv_CoDa_summary_epan = round(sum(colMeans(KLdiv_CoDa_epan)), 4) # 0.6998

## band_choice == "DPI"

# kernel = "gaussian"

KLdiv_CoDa_DPI = matrix(NA, 55, 2)
for(ik in 1:55)
{
    dat = cbind(true = test_density_DPI[,ik], forecast = den_fore_DPI[,ik])
    colnames(dat) = c("True", "Estimate")
    KLdiv_CoDa_DPI[ik,] = as.numeric(KLdiv(dat, eps=1e-16))[2:3]
    print(ik); rm(dat)
}

KLdiv_CoDa_summary_DPI = round(sum(colMeans(KLdiv_CoDa_DPI)), 4) # 0.8225

# kernel = "epan"

KLdiv_CoDa_epan_DPI = matrix(NA, 55, 2)
for(ik in 1:55)
{
    dat = cbind(true = test_density_epan_DPI[,ik], forecast = den_fore_epan_DPI[,ik])
    colnames(dat) = c("True", "Estimate")
    KLdiv_CoDa_epan_DPI[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
    print(ik); rm(dat)
}

KLdiv_CoDa_summary_epan_DPI = round(sum(colMeans(KLdiv_CoDa_epan_DPI)), 4) # 0.7495


##########################################
# Jensen-Shannon divergence (simple mean)
##########################################

## band_choice = "Silverman"

# kernel = "gaussian"

JSdiv_CoDa = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density[,ik], den_fore[,ik])
    M = rowMeans(CoDa_compar)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(den_fore[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_CoDa[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
}

JSdiv_CoDa_summary = round(sum(JSdiv_CoDa), 4) # 3.2215

# kernel = "epan"

JSdiv_CoDa_epan = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_epan[,ik], den_fore_epan[,ik])
    M = rowMeans(CoDa_compar)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(den_fore_epan[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_CoDa_epan[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
}

JSdiv_CoDa_summary_epan = round(sum(JSdiv_CoDa_epan), 4) # 2.3581

## band_choice = "DPI"

# kernel = "gaussian"

JSdiv_CoDa_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_DPI[,ik], den_fore_DPI[,ik])
    M = rowMeans(CoDa_compar)
    P_M = cbind(test_density_DPI[,ik], M)
    E_M = cbind(den_fore_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_CoDa_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
}

JSdiv_CoDa_summary_DPI = round(sum(JSdiv_CoDa_DPI), 4) # 3.835

# kernel = "epan"

JSdiv_CoDa_epan_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_epan_DPI[,ik], den_fore_epan_DPI[,ik])
    M = rowMeans(CoDa_compar)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(den_fore_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_CoDa_epan_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
}

JSdiv_CoDa_summary_epan_DPI = round(sum(JSdiv_CoDa_epan_DPI), 4) # 2.782

#############################################
# Jensen-Shannon divergence (geometric mean)
#############################################

## band_choice = "Silverman"

# kernel = "gaussian"

JSdiv_geo_CoDa = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density[,ik], den_fore[,ik])
    M = apply(CoDa_compar, 1, geometric.mean)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(den_fore[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_CoDa[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
}

JSdiv_geo_CoDa_summary = round(sum(JSdiv_geo_CoDa), 4) # 5.178

# kernel = "epan"

JSdiv_geo_CoDa_epan = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_epan[,ik], den_fore_epan[,ik])
    M = apply(CoDa_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(den_fore_epan[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_CoDa_epan[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(CoDa_compar)
}

JSdiv_geo_CoDa_epan_summary = round(sum(JSdiv_geo_CoDa_epan), 4) # 5.5992

## band_choice = "DPI"

# kernel = "gaussian"

JSdiv_geo_CoDa_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_DPI[,ik], den_fore_DPI[,ik])
    M = apply(CoDa_compar, 1, geometric.mean)
    P_M = cbind(test_density_DPI[,ik], M)
    E_M = cbind(den_fore_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_CoDa_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(CoDa_compar)
}

JSdiv_geo_CoDa_summary_DPI = round(sum(JSdiv_geo_CoDa_DPI), 4) # 6.4700

# kernel = "epan"

JSdiv_geo_CoDa_epan_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_epan_DPI[,ik], den_fore_epan_DPI[,ik])
    M = apply(CoDa_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(den_fore_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_CoDa_epan_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(CoDa_compar)
}

JSdiv_geo_CoDa_epan_summary_DPI = round(sum(JSdiv_geo_CoDa_epan_DPI), 4) # 6.0593


############
# L1 - norm
############

## band_choice = "Silverman"

# kernel = "gaussian"

L1_norm_CoDa = vector("numeric", 55)
for(ik in 1:55)
{
    L1_norm_CoDa[ik] = sum(abs(den_fore[,ik] - test_density[,ik]))
}

L1_norm_CoDa_summary = round(mean(L1_norm_CoDa), 4) # 951.5001

# kernel = "epan"

L1_norm_CoDa_epan = vector("numeric", 55)
for(ik in 1:55)
{
    L1_norm_CoDa_epan[ik] = sum(abs(den_fore_epan[,ik] - test_density_epan[,ik]))
}

L1_norm_CoDa_summary_epan = round(mean(L1_norm_CoDa_epan), 4) # 743.8618

## band_choice = "DPI"

# kernel = "gaussian"

L1_norm_CoDa_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    L1_norm_CoDa_DPI[ik] = sum(abs(den_fore_DPI[,ik] - test_density_DPI[,ik]))
}

L1_norm_CoDa_summary_DPI = round(mean(L1_norm_CoDa_DPI), 4) # 1055.207

# kernel = "epan"

L1_norm_CoDa_epan_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    L1_norm_CoDa_epan_DPI[ik] = sum(abs(den_fore_epan_DPI[,ik] - test_density_epan_DPI[,ik]))
}

L1_norm_CoDa_summary_epan_DPI = round(mean(L1_norm_CoDa_epan_DPI), 4) # 844.1048

############
# L2 - norm
############

## band_choice = "Silverman"

# kernel = "gaussian"

L2_norm_CoDa = vector("numeric", 55)
for(ik in 1:55)
{
    L2_norm_CoDa[ik] = sqrt(sum((test_density[,ik] - den_fore[,ik])^2))
}

L2_norm_CoDa_summary = round(mean(L2_norm_CoDa), 4) # 46.601

# kernel = "epan"

L2_norm_CoDa_epan = vector("numeric", 55)
for(ik in 1:55)
{
    L2_norm_CoDa_epan[ik] = sqrt(sum((test_density_epan[,ik] - den_fore_epan[,ik])^2))
}

L2_norm_CoDa_summary_epan = round(mean(L2_norm_CoDa_epan), 4) # 30.6467

## band_choice = "DPI"

# kernel = "gaussian"

L2_norm_CoDa_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    L2_norm_CoDa_DPI[ik] = sqrt(sum((test_density_DPI[,ik] - den_fore_DPI[,ik])^2))
}

L2_norm_CoDa_summary_DPI = round(mean(L2_norm_CoDa_DPI), 4) # 54.0202

# kernel = "epan"

L2_norm_CoDa_epan_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    L2_norm_CoDa_epan_DPI[ik] = sqrt(sum((test_density_epan_DPI[,ik] - den_fore_epan_DPI[,ik])^2))
}

L2_norm_CoDa_summary_epan_DPI = round(mean(L2_norm_CoDa_epan_DPI), 4) # 37.3269

##############
# Linf - norm
##############

## band_choice = "Silverman"

# kernel = "gaussian"

Linf_norm_CoDa = vector("numeric", 55)
for(ik in 1:55)
{
    Linf_norm_CoDa[ik] = max(abs(den_fore[,ik] - test_density[,ik]))
}

Linf_norm_CoDa_summary = round(mean(Linf_norm_CoDa), 4) # 3.6525

# kernel = "epan"

Linf_norm_CoDa_epan = vector("numeric", 55)
for(ik in 1:55)
{
    Linf_norm_CoDa_epan[ik] = max(abs(den_fore_epan[,ik] - test_density_epan[,ik]))    
}

Linf_norm_CoDa_summary_epan = round(mean(Linf_norm_CoDa_epan), 4) # 1.9194

## band_choice = "DPI"

# kernel = "gaussian"

Linf_norm_CoDa_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    Linf_norm_CoDa_DPI[ik] = max(abs(den_fore_DPI[,ik] - test_density_DPI[,ik]))
}

Linf_norm_CoDa_summary_DPI = round(mean(Linf_norm_CoDa_DPI), 4) # 4.4636

# kernel = "epan"

Linf_norm_CoDa_epan_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    Linf_norm_CoDa_epan_DPI[ik] = max(abs(den_fore_epan_DPI[,ik] - test_density_epan_DPI[,ik]))
}

Linf_norm_CoDa_summary_epan_DPI = round(mean(Linf_norm_CoDa_epan_DPI), 4) # 2.5715

# summary	

DJI_CoDa_summary = c(KLdiv_CoDa_summary,
                        JSdiv_CoDa_summary,
                        JSdiv_geo_CoDa_summary,
                        L1_norm_CoDa_summary,
                        L2_norm_CoDa_summary,
                        Linf_norm_CoDa_summary)

DJI_CoDa_summary_epan = c(KLdiv_CoDa_summary_epan,
                          JSdiv_CoDa_summary_epan,
                          JSdiv_geo_CoDa_epan_summary,
                          L1_norm_CoDa_summary_epan,
                          L2_norm_CoDa_summary_epan,
                          Linf_norm_CoDa_summary_epan)

DJI_CoDa_summary_DPI = c(KLdiv_CoDa_summary_DPI,
                         JSdiv_CoDa_summary_DPI,
                         JSdiv_geo_CoDa_summary_DPI,
                         L1_norm_CoDa_summary_DPI,
                         L2_norm_CoDa_summary_DPI,
                         Linf_norm_CoDa_summary_DPI)

DJI_CoDa_summary_epan_DPI = c(KLdiv_CoDa_summary_epan_DPI,
                              JSdiv_CoDa_summary_epan_DPI,
                              JSdiv_geo_CoDa_epan_summary_DPI,
                              L1_norm_CoDa_summary_epan_DPI,
                              L2_norm_CoDa_summary_epan_DPI,
                              Linf_norm_CoDa_summary_epan_DPI)

