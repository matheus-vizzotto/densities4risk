###########################
### CoDa method
### normalization = FALSE
###########################

## band_choice = "Silverman"

# kernel = "gaussian"

registerDoMC(8)
den_val_no_normalization = foreach(iw = 1:55) %dopar% DJI_forecast_den(data = DJI_ret, 
                                                                       normalization = FALSE, band_choice = "Silverman",
                                                                       ik = iw, kernel = "gaussian")

den_fore_no_normalization = matrix(NA, 5001, 55)
for(ik in 1:55)
{
    den_fore_no_normalization[,ik] = den_val_no_normalization[[ik]]
}
rm(den_val_no_normalization)

# kernel = "epanechnikov"

registerDoMC(8)
den_val_epan_no_normalization = foreach(iw = 1:55) %dopar% DJI_forecast_den(data = DJI_ret,
                                                                            normalization = FALSE, band_choice = "Silverman",
                                                                            ik = iw, kernel = "epanechnikov")

den_fore_epan_no_normalization = matrix(NA, 5001, 55)
for(ik in 1:55)
{
    den_fore_epan_no_normalization[,ik] = den_val_epan_no_normalization[[ik]]
}
rm(den_val_epan_no_normalization)

## band_choice = "DPI"

# kernel = "gaussian"

registerDoMC(8)
den_val_no_normalization_DPI = foreach(iw = 1:55) %dopar% DJI_forecast_den(data = DJI_ret, 
                                                                           normalization = FALSE, band_choice = "DPI",
                                                                           ik = iw, kernel = "gaussian")

den_fore_no_normalization_DPI = matrix(NA, 5001, 55)
for(ik in 1:55)
{
    den_fore_no_normalization_DPI[,ik] = den_val_no_normalization_DPI[[ik]]
}																															
rm(den_val_no_normalization_DPI)

# kernel = "epan"

registerDoMC(8)
den_val_epan_no_normalization_DPI = foreach(iw = 1:55) %dopar% DJI_forecast_den(data = DJI_ret,
                                                                                normalization = FALSE, band_choice = "DPI",
                                                                                ik = iw, kernel = "epanechnikov")

den_fore_epan_no_normalization_DPI = matrix(NA, 5001, 55)
for(ik in 1:55)
{
    den_fore_epan_no_normalization_DPI[,ik] = den_val_epan_no_normalization_DPI[[ik]]
}
rm(den_val_epan_no_normalization_DPI)


##############################
# Kullback-Leibler divergence
##############################

## band_choice = "Silverman"

# kernel = "gaussian"

KLdiv_CoDa_no_normalization = matrix(NA, 55, 2)
for(ik in 1:55)
{
    dat = cbind(true = test_density[,ik], forecast = den_fore_no_normalization[,ik])
    colnames(dat) = c("True", "Estimate")
    KLdiv_CoDa_no_normalization[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
    print(ik); rm(dat)
}

KLdiv_CoDa_no_normalization_summary = round(sum(colMeans(KLdiv_CoDa_no_normalization)), 4) # 0.6510

# kernel = "epan"

KLdiv_CoDa_epan_no_normalization = matrix(NA, 55, 2)
for(ik in 1:55)
{
    dat = cbind(true = test_density_epan[,ik], forecast = den_fore_epan_no_normalization[,ik])
    colnames(dat) = c("True", "Estimate")
    KLdiv_CoDa_epan_no_normalization[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
    print(ik); rm(dat)
}

KLdiv_CoDa_epan_no_normalization_summary = round(sum(colMeans(KLdiv_CoDa_epan_no_normalization)), 4) # 0.7004

## band_choice = "DPI"

# kernel = "gaussian"

KLdiv_CoDa_no_normalization_DPI = matrix(NA, 55, 2)
for(ik in 1:55)
{
    dat = cbind(true = test_density_DPI[,ik], forecast = den_fore_no_normalization_DPI[,ik])
    colnames(dat) = c("True", "Estimate")
    KLdiv_CoDa_no_normalization_DPI[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
    print(ik); rm(dat)
}

KLdiv_CoDa_no_normalization_summary_DPI = round(sum(colMeans(KLdiv_CoDa_no_normalization_DPI)), 4) # 0.8088

# kernel = "epan"

KLdiv_CoDa_epan_no_normalization_DPI = matrix(NA, 55, 2)
for(ik in 1:55)
{
    dat = cbind(true = test_density_epan_DPI[,ik], forecast = den_fore_epan_no_normalization_DPI[,ik])
    colnames(dat) = c("True", "Estimate")
    KLdiv_CoDa_epan_no_normalization_DPI[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
    print(ik); rm(dat)
}

KLdiv_CoDa_epan_no_normalization_summary_DPI = round(sum(colMeans(KLdiv_CoDa_epan_no_normalization_DPI)), 4) # 0.7476

#########################################
# Jensen-Shannon divergence (simple mean)
#########################################

## band_choice = "Silverman"

# kernel = "gaussian"

JSdiv_CoDa_no_normalization = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density[,ik], den_fore_no_normalization[,ik])
    M = rowMeans(CoDa_compar)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(den_fore_no_normalization[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_CoDa_no_normalization[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    print(ik); rm(CoDa_compar)
}

JSdiv_CoDa_no_normalization_summary = round(sum(JSdiv_CoDa_no_normalization), 4) # 3.1785

# kernel = "epan"

JSdiv_CoDa_epan_no_normalization = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_epan[,ik], den_fore_epan_no_normalization[,ik])
    M = rowMeans(CoDa_compar)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(den_fore_epan_no_normalization[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_CoDa_epan_no_normalization[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    print(ik); rm(CoDa_compar)
}

JSdiv_CoDa_epan_no_normalization_summary = round(sum(JSdiv_CoDa_epan_no_normalization), 4) # 2.3397

## band_choice = "DPI"

# kernel = "gaussian"

JSdiv_CoDa_no_normalization_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_DPI[,ik], den_fore_no_normalization_DPI[,ik])
    M = rowMeans(CoDa_compar)
    P_M = cbind(test_density_DPI[,ik], M)
    E_M = cbind(den_fore_no_normalization_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_CoDa_no_normalization_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    print(ik); rm(CoDa_compar)
}

JSdiv_CoDa_no_normalization_DPI_summary = round(sum(JSdiv_CoDa_no_normalization_DPI), 4) # 3.7884

# kernel = "epan"

JSdiv_CoDa_epan_no_normalization_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_epan_DPI[,ik], den_fore_epan_no_normalization_DPI[,ik])
    M = rowMeans(CoDa_compar)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(den_fore_epan_no_normalization_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_CoDa_epan_no_normalization_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    print(ik); rm(CoDa_compar)
}

JSdiv_CoDa_epan_no_normalization_DPI_summary = round(sum(JSdiv_CoDa_epan_no_normalization_DPI), 4) # 2.7679

############################################
# Jensen-Shannon divergence (geometric mean)
############################################

## band_choice = "Silverman"

# kernel = "gaussian"

JSdiv_geo_CoDa_no_normalization = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density[,ik], den_fore_no_normalization[,ik])
    M = apply(CoDa_compar, 1, geometric.mean)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(den_fore_no_normalization[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_CoDa_no_normalization[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    print(ik); rm(CoDa_compar)
}

JSdiv_geo_CoDa_no_normalization_summary = round(sum(JSdiv_geo_CoDa_no_normalization), 4) # 5.0572

# kernel = "epan"

JSdiv_geo_CoDa_epan_no_normalization = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_epan[,ik], den_fore_epan_no_normalization[,ik])
    M = apply(CoDa_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(den_fore_epan_no_normalization[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_CoDa_epan_no_normalization[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(CoDa_compar)
}

JSdiv_geo_CoDa_epan_no_normalization_summary = round(sum(JSdiv_geo_CoDa_epan_no_normalization), 4) # 5.5767

## band_choice = "DPI"

# kernel = "gaussian"

JSdiv_geo_CoDa_no_normalization_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_DPI[,ik], den_fore_no_normalization_DPI[,ik])
    M = apply(CoDa_compar, 1, geometric.mean)
    P_M = cbind(test_density_DPI[,ik], M)
    E_M = cbind(den_fore_no_normalization_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_CoDa_no_normalization_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(CoDa_compar)
}

JSdiv_geo_CoDa_no_normalization_DPI_summary = round(sum(JSdiv_geo_CoDa_no_normalization_DPI), 4) # 6.3581

# kernel = "epan"

JSdiv_geo_CoDa_epan_no_normalization_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    CoDa_compar = cbind(test_density_epan_DPI[,ik], den_fore_epan_no_normalization_DPI[,ik])
    M = apply(CoDa_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(den_fore_epan_no_normalization_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_CoDa_epan_no_normalization_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M); rm(CoDa_compar)
}

JSdiv_geo_CoDa_epan_no_normalization_DPI_summary = round(sum(JSdiv_geo_CoDa_epan_no_normalization_DPI), 4) # 6.0312

############
# L1 - norm
############

## band_choice = "Silverman"

# kernel = "gaussian"

L1_norm_CoDa_no_normalization = vector("numeric", 55)
for(ik in 1:55)
{
    L1_norm_CoDa_no_normalization[ik] = sum(abs(den_fore_no_normalization[,ik] - test_density[,ik]))
}

L1_norm_CoDa_no_normalization_summary = round(mean(L1_norm_CoDa_no_normalization), 4) # 943.6211

# kernel = "epan"

L1_norm_CoDa_epan_no_normalization = vector("numeric", 55)
for(ik in 1:55)
{
    L1_norm_CoDa_epan_no_normalization[ik] = sum(abs(den_fore_epan_no_normalization[,ik] - test_density_epan[,ik]))
}

L1_norm_CoDa_epan_no_normalization_summary = round(mean(L1_norm_CoDa_epan_no_normalization), 4) # 740.2241

## band_choice = "DPI"

# kernel = "gaussian"

L1_norm_CoDa_no_normalization_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    L1_norm_CoDa_no_normalization_DPI[ik] = sum(abs(den_fore_no_normalization_DPI[,ik] - test_density_DPI[,ik]))
}

L1_norm_CoDa_no_normalization_DPI_summary = round(mean(L1_norm_CoDa_no_normalization_DPI), 4) # 1049.785

# kernel = "epan"

L1_norm_CoDa_epan_no_normalization_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    L1_norm_CoDa_epan_no_normalization_DPI[ik] = sum(abs(den_fore_epan_no_normalization_DPI[,ik] - test_density_epan_DPI[,ik]))
}

L1_norm_CoDa_epan_no_normalization_DPI_summary = round(mean(L1_norm_CoDa_epan_no_normalization_DPI), 4) # 845.9398

############
# L2 - norm
############

## band_choice = "Silverman"

# kernel = "gaussian"

L2_norm_CoDa_no_normalization = vector("numeric", 55)
for(ik in 1:55)
{
    L2_norm_CoDa_no_normalization[ik] = sqrt(sum((test_density[,ik] - den_fore_no_normalization[,ik])^2))
}

L2_norm_CoDa_no_normalization_summary = round(mean(L2_norm_CoDa_no_normalization), 4) # 46.1948

# kernel = "epan"

L2_norm_CoDa_epan_no_normalization = vector("numeric", 55)
for(ik in 1:55)
{
    L2_norm_CoDa_epan_no_normalization[ik] = sqrt(sum((test_density_epan[,ik] - den_fore_epan_no_normalization[,ik])^2))
}

L2_norm_CoDa_epan_no_normalization_summary = round(mean(L2_norm_CoDa_epan_no_normalization), 4) # 30.4659

## band_choice = "DPI"

# kernel = "gaussian"

L2_norm_CoDa_no_normalization_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    L2_norm_CoDa_no_normalization_DPI[ik] = sqrt(sum((test_density_DPI[,ik] - den_fore_no_normalization_DPI[,ik])^2))
}

L2_norm_CoDa_no_normalization_DPI_summary = round(mean(L2_norm_CoDa_no_normalization_DPI), 4) # 53.7408

# kernel = "epan"

L2_norm_CoDa_epan_no_normalization_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    L2_norm_CoDa_epan_no_normalization_DPI[ik] = sqrt(sum((test_density_epan_DPI[,ik] - den_fore_epan_no_normalization_DPI[,ik])^2))
}

L2_norm_CoDa_epan_no_normalization_DPI_summary = round(mean(L2_norm_CoDa_epan_no_normalization_DPI), 4) # 37.4132

##############
# Linf - norm
##############

## band_choice = "Silverman"

# kernel = "gaussian"

Linf_norm_CoDa_no_normalization = vector("numeric", 55)
for(ik in 1:55)
{
    Linf_norm_CoDa_no_normalization[ik] = max(abs(den_fore_no_normalization[,ik] - test_density[,ik]))
}

Linf_norm_CoDa_no_normalization_summary = round(mean(Linf_norm_CoDa_no_normalization), 4) # 3.6227

# kernel = "epan"

Linf_norm_CoDa_epan_no_normalization = vector("numeric", 55)
for(ik in 1:55)
{
    Linf_norm_CoDa_epan_no_normalization[ik] = max(abs(den_fore_epan_no_normalization[,ik] - test_density_epan[,ik]))
}

Linf_norm_CoDa_epan_no_normalization_summary = round(mean(Linf_norm_CoDa_epan_no_normalization), 4) # 1.9031

## band_choice = "DPI"

# kernel = "gaussian"

Linf_norm_CoDa_no_normalization_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    Linf_norm_CoDa_no_normalization_DPI[ik] = max(abs(den_fore_no_normalization_DPI[,ik] - test_density_DPI[,ik]))
}

Linf_norm_CoDa_no_normalization_DPI_summary = round(mean(Linf_norm_CoDa_no_normalization_DPI), 4) # 4.4341

# kernel = "epan"

Linf_norm_CoDa_epan_no_normalization_DPI = vector("numeric", 55)
for(ik in 1:55)
{
    Linf_norm_CoDa_epan_no_normalization_DPI[ik] = max(abs(den_fore_epan_no_normalization_DPI[,ik] - test_density_epan_DPI[,ik]))
}

Linf_norm_CoDa_epan_no_normalization_DPI_summary = round(mean(Linf_norm_CoDa_epan_no_normalization_DPI), 4) # 2.5866

# summary

DJI_CoDa_no_normalization_summary = c(KLdiv_CoDa_no_normalization_summary,
                                      JSdiv_CoDa_no_normalization_summary,
                                      JSdiv_geo_CoDa_no_normalization_summary,
                                      L1_norm_CoDa_no_normalization_summary,
                                      L2_norm_CoDa_no_normalization_summary,
                                      Linf_norm_CoDa_no_normalization_summary)

DJI_CoDa_epan_no_normalization_summary = c(KLdiv_CoDa_epan_no_normalization_summary,
                                           JSdiv_CoDa_epan_no_normalization_summary,
                                           JSdiv_geo_CoDa_epan_no_normalization_summary,
                                           L1_norm_CoDa_epan_no_normalization_summary,
                                           L2_norm_CoDa_epan_no_normalization_summary,
                                           Linf_norm_CoDa_epan_no_normalization_summary)

DJI_CoDa_no_normalization_DPI_summary = c(KLdiv_CoDa_no_normalization_summary_DPI,
                                          JSdiv_CoDa_no_normalization_DPI_summary,
                                          JSdiv_geo_CoDa_no_normalization_DPI_summary,
                                          L1_norm_CoDa_no_normalization_DPI_summary,
                                          L2_norm_CoDa_no_normalization_DPI_summary,
                                          Linf_norm_CoDa_no_normalization_DPI_summary)

DJI_CoDa_epan_no_normalization_DPI_summary = c(KLdiv_CoDa_epan_no_normalization_summary_DPI,
                                            JSdiv_CoDa_epan_no_normalization_DPI_summary,
                                            JSdiv_geo_CoDa_epan_no_normalization_DPI_summary,
                                            L1_norm_CoDa_epan_no_normalization_DPI_summary,
                                            L2_norm_CoDa_epan_no_normalization_DPI_summary,
                                            Linf_norm_CoDa_epan_no_normalization_DPI_summary)

