######################################
# Forecasts (univariate FTS approach)
######################################

# fmethod = "ets"

Alex_ets_DJI_ret_fore_den_DPI = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_ets_DJI_ret_fore_den_DPI[,ik] = Alex(data = DJI_ret[1:(109+ik),], forecasting_method = "uni", 
                                              band_choice = "DPI", fmethod = "ets", kernel = "gaussian")$dens_fore
    print(109+ik)
}

Alex_ets_DJI_ret_fore_den_epan_DPI = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_ets_DJI_ret_fore_den_epan_DPI[,ik] = Alex(data = DJI_ret[1:(109+ik),], forecasting_method = "uni", 
                                                   band_choice = "DPI", fmethod = "ets", kernel = "epanechnikov")$dens_fore
    print(109+ik)
}

# fmethod = "arima"

Alex_arima_DJI_ret_fore_den_DPI = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_arima_DJI_ret_fore_den_DPI[,ik] = Alex(data = DJI_ret[1:(109+ik),], forecasting_method = "uni", 
                                                band_choice = "DPI", fmethod = "arima", kernel = "gaussian")$dens_fore
    print(109+ik)
}

Alex_arima_DJI_ret_fore_den_epan_DPI = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_arima_DJI_ret_fore_den_epan_DPI[,ik] = Alex(data = DJI_ret[1:(109+ik),], forecasting_method = "uni", 
                                                     band_choice = "DPI", fmethod = "arima", kernel = "epanechnikov")$dens_fore
    print(109+ik)
}

###########################
# Forecasts (BYZ approach)
###########################

Alex_BYZ_DJI_ret_fore_den_DPI = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_BYZ_DJI_ret_fore_den_DPI[,ik] = Alex(data = DJI_ret[1:(109+ik),], forecasting_method = "BYZ", 
                                              band_choice = "DPI", kernel = "gaussian")$dens_fore
    print(ik)
}

Alex_BYZ_DJI_ret_fore_den_epan_DPI = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Alex_BYZ_DJI_ret_fore_den_epan_DPI[,ik] = Alex(data = DJI_ret[1:(109+ik),], forecasting_method = "BYZ", 
                                                   band_choice = "DPI", kernel = "epanechnikov")$dens_fore
    print(109+ik)
}

##############################
# Kullback-Leibler divergence
##############################

# kernel = "gaussian"

KLdiv_Alex_ets_DPI = KLdiv_Alex_arima_DPI = KLdiv_Alex_BYZ_DPI = matrix(NA, n_test, 2)
for(ik in 1:n_test)
{
    Alex_ets_compar = cbind(test_density_DPI[,ik], Alex_ets_DJI_ret_fore_den_DPI[,ik])
    Alex_arima_compar = cbind(test_density_DPI[,ik], Alex_arima_DJI_ret_fore_den_DPI[,ik])
    Alex_BYZ_compar = cbind(test_density_DPI[,ik], Alex_BYZ_DJI_ret_fore_den_DPI[,ik])
    colnames(Alex_ets_compar) = colnames(Alex_arima_compar) = colnames(Alex_BYZ_compar) = c("True", "Estimate")
    
    KLdiv_Alex_ets_DPI[ik,] = as.numeric(KLdiv(Alex_ets_compar, eps = 1e-16))[2:3]
    KLdiv_Alex_arima_DPI[ik,] = as.numeric(KLdiv(Alex_arima_compar, eps = 1e-16))[2:3]
    KLdiv_Alex_BYZ_DPI[ik,] = as.numeric(KLdiv(Alex_BYZ_compar, eps = 1e-16))[2:3]
    print(ik)
    rm(Alex_ets_compar); rm(Alex_arima_compar); rm(Alex_BYZ_compar)
}

KLdiv_Alex_ets_DPI_summary = round(sum(colMeans(KLdiv_Alex_ets_DPI)), 4) # 1.3165
KLdiv_Alex_arima_DPI_summary = round(sum(colMeans(KLdiv_Alex_arima_DPI)), 4) # 1.2569
KLdiv_Alex_BYZ_DPI_summary = round(sum(colMeans(KLdiv_Alex_BYZ_DPI)), 4) # 1.3312

# kernel = "epan"

KLdiv_Alex_ets_epan_DPI = KLdiv_Alex_arima_epan_DPI = KLdiv_Alex_BYZ_epan_DPI = matrix(NA, n_test, 2)
for(ik in 1:n_test)
{
    Alex_ets_compar = cbind(test_density_epan_DPI[,ik], Alex_ets_DJI_ret_fore_den_epan_DPI[,ik])
    Alex_arima_compar = cbind(test_density_epan_DPI[,ik], Alex_arima_DJI_ret_fore_den_epan_DPI[,ik])
    Alex_BYZ_compar = cbind(test_density_epan_DPI[,ik], Alex_BYZ_DJI_ret_fore_den_epan_DPI[,ik])
    colnames(Alex_ets_compar) = colnames(Alex_arima_compar) = colnames(Alex_BYZ_compar) = c("True", "Estimate")
    
    KLdiv_Alex_ets_epan_DPI[ik,] = as.numeric(KLdiv(Alex_ets_compar, eps = 1e-16))[2:3]
    KLdiv_Alex_arima_epan_DPI[ik,] = as.numeric(KLdiv(Alex_arima_compar, eps = 1e-16))[2:3]
    KLdiv_Alex_BYZ_epan_DPI[ik,] = as.numeric(KLdiv(Alex_BYZ_compar, eps = 1e-16))[2:3]
    print(ik)
    rm(Alex_ets_compar); rm(Alex_arima_compar); rm(Alex_BYZ_compar)
}

KLdiv_Alex_ets_epan_DPI_summary = round(sum(colMeans(KLdiv_Alex_ets_epan_DPI)), 4) # 1.1629
KLdiv_Alex_arima_epan_DPI_summary = round(sum(colMeans(KLdiv_Alex_arima_epan_DPI)), 4) # 1.1197
KLdiv_Alex_BYZ_epan_DPI_summary = round(sum(colMeans(KLdiv_Alex_BYZ_epan_DPI)), 4) # 1.1485

############################
# Jensen-Shannon divergence
############################

# kernel = "gaussian"

JSdiv_Alex_ets_DPI = JSdiv_Alex_arima_DPI = JSdiv_Alex_BYZ_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Alex_ets_compar = cbind(test_density_DPI[,ik], Alex_ets_DJI_ret_fore_den_DPI[,ik])
    M = rowMeans(Alex_ets_compar)
    P_M = cbind(test_density_DPI[,ik], M)
    E_M = cbind(Alex_ets_DJI_ret_fore_den_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_ets_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_ets_compar)
    
    Alex_arima_compar = cbind(test_density_DPI[,ik], Alex_arima_DJI_ret_fore_den_DPI[,ik])
    M = rowMeans(Alex_arima_compar)
    P_M = cbind(test_density_DPI[,ik], M)
    E_M = cbind(Alex_arima_DJI_ret_fore_den_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_arima_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_arima_compar)
    
    Alex_BYZ_compar = cbind(test_density_DPI[,ik], Alex_BYZ_DJI_ret_fore_den_DPI[,ik])
    M = rowMeans(Alex_BYZ_compar)
    P_M = cbind(test_density_DPI[,ik], M)
    E_M = cbind(Alex_BYZ_DJI_ret_fore_den_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_BYZ_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_BYZ_compar)
    print(ik)
}

JSdiv_Alex_ets_DPI_summary = round(sum(JSdiv_Alex_ets_DPI), 4) # 3.5591
JSdiv_Alex_arima_DPI_summary = round(sum(JSdiv_Alex_arima_DPI), 4) # 3.4896
JSdiv_Alex_BYZ_DPI_summary = round(sum(JSdiv_Alex_BYZ_DPI), 4) # 3.4860

# kernel = "epan"

JSdiv_Alex_ets_epan_DPI = JSdiv_Alex_arima_epan_DPI = JSdiv_Alex_BYZ_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Alex_ets_compar = cbind(test_density_epan_DPI[,ik], Alex_ets_DJI_ret_fore_den_epan_DPI[,ik])
    M = rowMeans(Alex_ets_compar)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(Alex_ets_DJI_ret_fore_den_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_ets_epan_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_ets_compar)
    
    Alex_arima_compar = cbind(test_density_epan_DPI[,ik], Alex_arima_DJI_ret_fore_den_epan_DPI[,ik])
    M = rowMeans(Alex_arima_compar)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(Alex_arima_DJI_ret_fore_den_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_arima_epan_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_arima_compar)
    
    Alex_BYZ_compar = cbind(test_density_epan_DPI[,ik], Alex_BYZ_DJI_ret_fore_den_epan_DPI[,ik])
    M = rowMeans(Alex_BYZ_compar)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(Alex_BYZ_DJI_ret_fore_den_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Alex_BYZ_epan_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_BYZ_compar)
    print(ik)
}

JSdiv_Alex_ets_epan_DPI_summary = round(sum(JSdiv_Alex_ets_epan_DPI), 4) # 2.1894
JSdiv_Alex_arima_epan_DPI_summary = round(sum(JSdiv_Alex_arima_epan_DPI), 4) # 2.1370
JSdiv_Alex_BYZ_epan_DPI_summary = round(sum(JSdiv_Alex_BYZ_epan_DPI), 4) # 2.1362

##################################
# Jensen-Shannon divergence (geo)
##################################

# kernel = "gaussian"

JSdiv_geo_Alex_ets_DPI = JSdiv_geo_Alex_arima_DPI = JSdiv_geo_Alex_BYZ_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Alex_ets_compar = cbind(test_density_DPI[,ik], Alex_ets_DJI_ret_fore_den_epan_DPI[,ik])
    M = apply(Alex_ets_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(Alex_ets_DJI_ret_fore_den_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_ets_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_ets_compar)
    
    Alex_arima_compar = cbind(test_density_epan_DPI[,ik], Alex_arima_DJI_ret_fore_den_epan_DPI[,ik])
    M = apply(Alex_arima_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(Alex_arima_DJI_ret_fore_den_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_arima_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_arima_compar)
    
    Alex_BYZ_compar = cbind(test_density_epan_DPI[,ik], Alex_BYZ_DJI_ret_fore_den_epan_DPI[,ik])
    M = apply(Alex_BYZ_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(Alex_BYZ_DJI_ret_fore_den_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_BYZ_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_BYZ_compar)
    print(ik)
}

JSdiv_geo_Alex_ets_DPI_summary = round(sum(JSdiv_geo_Alex_ets_DPI), 4) # 10.9937
JSdiv_geo_Alex_arima_DPI_summary = round(sum(JSdiv_geo_Alex_arima_DPI), 4) # 3.4896
JSdiv_geo_Alex_BYZ_DPI_summary = round(sum(JSdiv_geo_Alex_BYZ_DPI), 4) # 3.4860

# kernel = "epan"

JSdiv_geo_Alex_ets_epan_DPI = JSdiv_geo_Alex_arima_epan_DPI = JSdiv_geo_Alex_BYZ_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Alex_ets_compar = cbind(test_density_epan_DPI[,ik], Alex_ets_DJI_ret_fore_den_epan_DPI[,ik])
    M = apply(Alex_ets_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(Alex_ets_DJI_ret_fore_den_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_ets_epan_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_ets_compar)
    
    Alex_arima_compar = cbind(test_density_epan_DPI[,ik], Alex_arima_DJI_ret_fore_den_epan_DPI[,ik])
    M = apply(Alex_arima_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(Alex_arima_DJI_ret_fore_den_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_arima_epan_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_arima_compar)
    
    Alex_BYZ_compar = cbind(test_density_epan_DPI[,ik], Alex_BYZ_DJI_ret_fore_den_epan_DPI[,ik])
    M = apply(Alex_BYZ_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(Alex_BYZ_DJI_ret_fore_den_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Alex_BYZ_epan_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(Alex_BYZ_compar)
    print(ik)
}

JSdiv_geo_Alex_ets_epan_DPI_summary = round(sum(JSdiv_geo_Alex_ets_epan_DPI), 4) # 6.9713
JSdiv_geo_Alex_arima_epan_DPI_summary = round(sum(JSdiv_geo_Alex_arima_epan_DPI), 4) # 6.7053
JSdiv_geo_Alex_BYZ_epan_DPI_summary = round(sum(JSdiv_geo_Alex_BYZ_epan_DPI), 4) # 6.8762

############
# L1 - norm
############

# kernel = "gaussian"

L1_norm_Alex_ets_DPI = L1_norm_Alex_arima_DPI = L1_norm_Alex_BYZ_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_Alex_ets_DPI[ik] = sum(abs(Alex_ets_DJI_ret_fore_den_DPI[,ik] - test_density_DPI[,ik]))
    L1_norm_Alex_arima_DPI[ik] = sum(abs(Alex_arima_DJI_ret_fore_den_DPI[,ik] - test_density_DPI[,ik]))
    L1_norm_Alex_BYZ_DPI[ik] = sum(abs(Alex_BYZ_DJI_ret_fore_den_DPI[,ik] - test_density_DPI[,ik]))
}

L1_norm_Alex_ets_DPI_summary = round(mean(L1_norm_Alex_ets_DPI), 4) # 1060.639
L1_norm_Alex_arima_DPI_summary = round(mean(L1_norm_Alex_arima_DPI), 4) # 1043.122
L1_norm_Alex_BYZ_DPI_summary = round(mean(L1_norm_Alex_BYZ_DPI), 4) # 1040.163

# kernel = "epan"

L1_norm_Alex_ets_epan_DPI = L1_norm_Alex_arima_epan_DPI = L1_norm_Alex_BYZ_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_Alex_ets_epan_DPI[ik] = sum(abs(Alex_ets_DJI_ret_fore_den_epan_DPI[,ik] - test_density_epan_DPI[,ik]))
    L1_norm_Alex_arima_epan_DPI[ik] = sum(abs(Alex_arima_DJI_ret_fore_den_epan_DPI[,ik] - test_density_epan_DPI[,ik]))
    L1_norm_Alex_BYZ_epan_DPI[ik] = sum(abs(Alex_BYZ_DJI_ret_fore_den_epan_DPI[,ik] - test_density_epan_DPI[,ik]))
}

L1_norm_Alex_ets_epan_DPI_summary = round(mean(L1_norm_Alex_ets_epan_DPI), 4) # 772.90
L1_norm_Alex_arima_epan_DPI_summary = round(mean(L1_norm_Alex_arima_epan_DPI), 4) # 757.2343
L1_norm_Alex_BYZ_epan_DPI_summary = round(mean(L1_norm_Alex_BYZ_epan_DPI), 4) # 754.2996

############
# L2 - norm
############

# kernel = "gaussian"

L2_norm_Alex_ets_DPI = L2_norm_Alex_arima_DPI = L2_norm_Alex_BYZ_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_Alex_ets_DPI[ik] = sqrt(sum((test_density_DPI[,ik] - Alex_ets_DJI_ret_fore_den_DPI[,ik])^2))
    L2_norm_Alex_arima_DPI[ik] = sqrt(sum((test_density_DPI[,ik] - Alex_arima_DJI_ret_fore_den_DPI[,ik])^2))
    L2_norm_Alex_BYZ_DPI[ik] = sqrt(sum((test_density_DPI[,ik] - Alex_BYZ_DJI_ret_fore_den_DPI[,ik])^2))
}

L2_norm_Alex_ets_DPI_summary = round(mean(L2_norm_Alex_ets_DPI), 4) # 53.2799
L2_norm_Alex_arima_DPI_summary = round(mean(L2_norm_Alex_arima_DPI), 4) # 52.3588
L2_norm_Alex_BYZ_DPI_summary = round(mean(L2_norm_Alex_BYZ_DPI), 4) # 52.1259

# kernel = "epan"

L2_norm_Alex_ets_epan_DPI = L2_norm_Alex_arima_epan_DPI = L2_norm_Alex_BYZ_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_Alex_ets_epan_DPI[ik] = sqrt(sum((test_density_epan_DPI[,ik] - Alex_ets_DJI_ret_fore_den_epan_DPI[,ik])^2))
    L2_norm_Alex_arima_epan_DPI[ik] = sqrt(sum((test_density_epan_DPI[,ik] - Alex_arima_DJI_ret_fore_den_epan_DPI[,ik])^2))
    L2_norm_Alex_BYZ_epan_DPI[ik] = sqrt(sum((test_density_epan_DPI[,ik] - Alex_BYZ_DJI_ret_fore_den_epan_DPI[,ik])^2))
}

L2_norm_Alex_ets_epan_DPI_summary = round(mean(L2_norm_Alex_ets_epan_DPI), 4) # 33.3286
L2_norm_Alex_arima_epan_DPI_summary = round(mean(L2_norm_Alex_arima_epan_DPI), 4) # 32.6422
L2_norm_Alex_BYZ_epan_DPI_summary = round(mean(L2_norm_Alex_BYZ_epan_DPI), 4) # 32.4956

##############
# Linf - norm
##############

# kernel = "gaussian"

Linf_norm_Alex_ets_DPI = Linf_norm_Alex_arima_DPI = Linf_norm_Alex_BYZ_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_Alex_ets_DPI[ik] = max(abs(Alex_ets_DJI_ret_fore_den_DPI[,ik] - test_density_DPI[,ik]))
    Linf_norm_Alex_arima_DPI[ik] = max(abs(Alex_arima_DJI_ret_fore_den_DPI[,ik] - test_density_DPI[,ik]))
    Linf_norm_Alex_BYZ_DPI[ik] = max(abs(Alex_BYZ_DJI_ret_fore_den_DPI[,ik] - test_density_DPI[,ik]))
}

Linf_norm_Alex_ets_DPI_summary = round(mean(Linf_norm_Alex_ets_DPI), 4) # 4.4753
Linf_norm_Alex_arima_DPI_summary = round(mean(Linf_norm_Alex_arima_DPI), 4) # 4.4017
Linf_norm_Alex_BYZ_DPI_summary = round(mean(Linf_norm_Alex_BYZ_DPI), 4) # 4.3670

# kernel = "epan"

Linf_norm_Alex_ets_epan_DPI = Linf_norm_Alex_arima_epan_DPI = Linf_norm_Alex_BYZ_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_Alex_ets_epan_DPI[ik] = max(abs(Alex_ets_DJI_ret_fore_den_epan_DPI[,ik] - test_density_epan_DPI[,ik]))
    Linf_norm_Alex_arima_epan_DPI[ik] = max(abs(Alex_arima_DJI_ret_fore_den_epan_DPI[,ik] - test_density_epan_DPI[,ik]))
    Linf_norm_Alex_BYZ_epan_DPI[ik] = max(abs(Alex_BYZ_DJI_ret_fore_den_epan_DPI[,ik] - test_density_epan_DPI[,ik]))
}

Linf_norm_Alex_ets_epan_DPI_summary = round(mean(Linf_norm_Alex_ets_epan_DPI), 4) # 2.2189
Linf_norm_Alex_arima_epan_DPI_summary = round(mean(Linf_norm_Alex_arima_epan_DPI), 4) # 2.1657
Linf_norm_Alex_BYZ_epan_DPI_summary = round(mean(Linf_norm_Alex_BYZ_epan_DPI), 4) # 2.1401

##########
# summary
##########

DJI_LQD_arima_DPI_summary = c(KLdiv_Alex_arima_DPI_summary,
                                JSdiv_Alex_arima_DPI_summary,
                                JSdiv_geo_Alex_arima_DPI_summary,
                                L1_norm_Alex_arima_DPI_summary,
                                L2_norm_Alex_arima_DPI_summary,
                                Linf_norm_Alex_arima_DPI_summary)

DJI_LQD_arima_DPI_epan_summary = c(KLdiv_Alex_arima_epan_DPI_summary,
                                    JSdiv_Alex_arima_epan_DPI_summary,
                                    JSdiv_geo_Alex_arima_epan_DPI_summary,
                                    L1_norm_Alex_arima_epan_DPI_summary,
                                    L2_norm_Alex_arima_epan_DPI_summary,
                                    Linf_norm_Alex_arima_epan_DPI_summary)

DJI_LQD_ets_DPI_summary = c(KLdiv_Alex_ets_DPI_summary,
                            JSdiv_Alex_ets_DPI_summary,
                            JSdiv_geo_Alex_ets_DPI_summary,
                            L1_norm_Alex_ets_DPI_summary,
                            L2_norm_Alex_ets_DPI_summary,
                            Linf_norm_Alex_ets_DPI_summary)

DJI_LQD_ets_DPI_epan_summary = c(KLdiv_Alex_ets_epan_DPI_summary,
                                JSdiv_Alex_ets_epan_DPI_summary,
                                JSdiv_geo_Alex_ets_epan_DPI_summary,
                                L1_norm_Alex_ets_epan_DPI_summary,
                                L2_norm_Alex_ets_epan_DPI_summary,
                                Linf_norm_Alex_ets_epan_DPI_summary)

DJI_LQD_BYZ_DPI_summary = c(KLdiv_Alex_BYZ_DPI_summary,
                            JSdiv_Alex_BYZ_DPI_summary,
                            JSdiv_geo_Alex_BYZ_DPI_summary,
                            L1_norm_Alex_BYZ_DPI_summary,
                            L2_norm_Alex_BYZ_DPI_summary,
                            Linf_norm_Alex_BYZ_DPI_summary)

DJI_LQD_BYZ_DPI_epan_summary = c(KLdiv_Alex_BYZ_epan_DPI_summary,
                                JSdiv_Alex_BYZ_epan_DPI_summary,
                                JSdiv_geo_Alex_BYZ_epan_DPI_summary,
                                L1_norm_Alex_BYZ_epan_DPI_summary,
                                L2_norm_Alex_BYZ_epan_DPI_summary,
                                Linf_norm_Alex_BYZ_epan_DPI_summary)

