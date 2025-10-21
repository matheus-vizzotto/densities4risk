###########
# Test set
###########

test_density_DPI = ken_density(data = DJI_ret, band_choice = "DPI", kernel = "gaussian")$Y[,111:165]

test_density_epan_DPI = ken_density(data = DJI_ret, band_choice = "DPI", kernel = "epanechnikov")$Y[,111:165]

###############################
# Horta-Ziegelman (auto select)
###############################

# kernel = "gaussian"

DJI_ret_fore_den_DPI = matrix(NA, 5001, n_test)
DJI_ret_fore_d0_DPI = vector("numeric", n_test)
for(ik in 1:55)
{
    dum = Horta_Ziegelman(data = DJI_ret[1:(109+ik),], band_choice = "DPI", kernel = "gaussian")
    DJI_ret_fore_den_DPI[,ik] = dum$Yhat.fix_den
    DJI_ret_fore_d0_DPI[ik] = dum$selected_d0
    print(109+ik); rm(dum); rm(ik)
}

# kernel = "epan"

DJI_ret_fore_den_epan_DPI = matrix(NA, 5001, n_test)
DJI_ret_fore_d0_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    dum = Horta_Ziegelman(data = DJI_ret[1:(109+ik),], band_choice = "DPI", kernel = "epanechnikov")
    DJI_ret_fore_den_epan_DPI[,ik] = dum$Yhat.fix_den
    DJI_ret_fore_d0_epan_DPI[ik] = dum$selected_d0
    print(109+ik); rm(dum); rm(ik)
}

#####################
# Density evaluation
#####################

## Kullback-Leibler divergence

# kernel = "gaussian"

KLdiv_Horta_Ziegelman_DPI = matrix(NA, n_test, 2)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density_DPI[,ik], DJI_ret_fore_den_DPI[,ik])
    colnames(Horta_Ziegelman_compar) = c("True", "Estimate")
    KLdiv_Horta_Ziegelman_DPI[ik,] = as.numeric(KLdiv(Horta_Ziegelman_compar, eps = 1e-16))[2:3]
    print(ik); rm(Horta_Ziegelman_compar)
}

round(colMeans(KLdiv_Horta_Ziegelman_DPI), 4) # 1.1973 0.2778
KLdiv_Horta_Ziegelman_DPI_summary = round(sum(colMeans(KLdiv_Horta_Ziegelman_DPI)), 4) # 1.4751

# kernel = "epan"

KLdiv_Horta_Ziegelman_epan_DPI = matrix(NA, n_test, 2)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density_epan_DPI[,ik], DJI_ret_fore_den_epan_DPI[,ik])
    colnames(Horta_Ziegelman_compar) = c("True", "Estimate")
    KLdiv_Horta_Ziegelman_epan_DPI[ik,] = as.numeric(KLdiv(Horta_Ziegelman_compar, eps = 1e-16))[2:3]
    print(ik); rm(Horta_Ziegelman_compar)
}

KLdiv_Horta_Ziegelman_epan_DPI_summary = round(sum(colMeans(KLdiv_Horta_Ziegelman_epan_DPI)), 4) # 1.7475

## Jensen-Shannon divergence

# kernel = "gaussian"

JSdiv_Horta_Ziegelman_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density_DPI[,ik], DJI_ret_fore_den_DPI[,ik])
    M = rowMeans(Horta_Ziegelman_compar)
    P_M = cbind(test_density_DPI[,ik], M)
    E_M = cbind(DJI_ret_fore_den_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Horta_Ziegelman_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    print(ik); rm(Horta_Ziegelman_compar)
}

JSdiv_Horta_Ziegelman_DPI_summary = round(sum(JSdiv_Horta_Ziegelman_DPI), 4) # 4.0955

# kernel = "epan"

JSdiv_Horta_Ziegelman_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density_epan_DPI[,ik], DJI_ret_fore_den_epan_DPI[,ik])
    M = rowMeans(Horta_Ziegelman_compar)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(DJI_ret_fore_den_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_Horta_Ziegelman_epan_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    print(ik); rm(Horta_Ziegelman_compar)
}

JSdiv_Horta_Ziegelman_epan_DPI_summary = round(sum(JSdiv_Horta_Ziegelman_epan_DPI), 4) # 2.5834

## Jensen-Shannon divergence (geometric mean)

# kernel = "gaussian"

JSdiv_geo_Horta_Ziegelman_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density_DPI[,ik], DJI_ret_fore_den_DPI[,ik])
    M = apply(Horta_Ziegelman_compar, 1, geometric.mean)
    P_M = cbind(test_density_DPI[,ik], M)
    E_M = cbind(DJI_ret_fore_den_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Horta_Ziegelman_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    print(ik); rm(Horta_Ziegelman_compar)
}

JSdiv_geo_Horta_Ziegelman_DPI_summary = round(sum(JSdiv_geo_Horta_Ziegelman_DPI), 4) # 10.6577

# kernel = "epan"

JSdiv_geo_Horta_Ziegelman_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Horta_Ziegelman_compar = cbind(test_density_epan_DPI[,ik], DJI_ret_fore_den_epan_DPI[,ik])
    M = apply(Horta_Ziegelman_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(DJI_ret_fore_den_epan_DPI[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    JSdiv_geo_Horta_Ziegelman_epan_DPI[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    print(ik); rm(Horta_Ziegelman_compar)
}

JSdiv_geo_Horta_Ziegelman_epan_DPI_summary = round(sum(JSdiv_geo_Horta_Ziegelman_epan_DPI), 4) # 9.8658

## L1-norm

# kernel = "gaussian"

L1_norm_Horta_Ziegelman_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_Horta_Ziegelman_DPI[ik] = sum(abs(test_density_DPI[,ik] - DJI_ret_fore_den_DPI[,ik]))
}

L1_norm_Horta_Ziegelman_DPI_summary = round(mean(L1_norm_Horta_Ziegelman_DPI), 4) # 1144.323

# kernel = "epan"

L1_norm_Horta_Ziegelman_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_Horta_Ziegelman_epan_DPI[ik] = sum(abs(test_density_epan_DPI[,ik] - DJI_ret_fore_den_epan_DPI[,ik]))
}

L1_norm_Horta_Ziegelman_epan_DPI_summary = round(mean(L1_norm_Horta_Ziegelman_epan_DPI), 4) # 835.1553

## L2-norm

# kernel = "gaussian"

L2_norm_Horta_Ziegelman_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_Horta_Ziegelman_DPI[ik] = sqrt(sum((test_density_DPI[,ik] - DJI_ret_fore_den_DPI[,ik])^2))
}

L2_norm_Horta_Ziegelman_DPI_summary = round(mean(L2_norm_Horta_Ziegelman_DPI), 4) # 53.7560

# kernel = "epan"

L2_norm_Horta_Ziegelman_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_Horta_Ziegelman_epan_DPI[ik] = sqrt(sum((test_density_epan_DPI[,ik] - DJI_ret_fore_den_epan_DPI[,ik])^2))
}

L2_norm_Horta_Ziegelman_epan_DPI_summary = round(mean(L2_norm_Horta_Ziegelman_epan_DPI), 4) # 33.7902

## Linf-norm

# kernel = "gaussian"

Linf_norm_Horta_Ziegelman_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_Horta_Ziegelman_DPI[ik] = max(abs(test_density_DPI[,ik] - DJI_ret_fore_den_DPI[,ik]))
}

Linf_norm_Horta_Ziegelman_DPI_summary = round(mean(Linf_norm_Horta_Ziegelman_DPI), 4) # 4.8392

# kernel = "epan"

Linf_norm_Horta_Ziegelman_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_Horta_Ziegelman_epan_DPI[ik] = max(abs(test_density_epan_DPI[,ik] - DJI_ret_fore_den_epan_DPI[,ik]))
}

Linf_norm_Horta_Ziegelman_epan_DPI_summary = round(mean(Linf_norm_Horta_Ziegelman_epan_DPI), 4) # 2.2728

##########
# summary
##########

DJI_Horta_Ziegelman_DPI_summary = c(KLdiv_Horta_Ziegelman_DPI_summary,
                                    JSdiv_Horta_Ziegelman_DPI_summary,
                                    JSdiv_geo_Horta_Ziegelman_DPI_summary,
                                    L1_norm_Horta_Ziegelman_DPI_summary,
                                    L2_norm_Horta_Ziegelman_DPI_summary,
                                    Linf_norm_Horta_Ziegelman_DPI_summary)

DJI_Horta_Ziegelman_epan_DPI_summary = c(KLdiv_Horta_Ziegelman_epan_DPI_summary,
                                         JSdiv_Horta_Ziegelman_epan_DPI_summary,
                                         JSdiv_geo_Horta_Ziegelman_epan_DPI_summary,
                                         L1_norm_Horta_Ziegelman_epan_DPI_summary,
                                         L2_norm_Horta_Ziegelman_epan_DPI_summary,
                                         Linf_norm_Horta_Ziegelman_epan_DPI_summary)

