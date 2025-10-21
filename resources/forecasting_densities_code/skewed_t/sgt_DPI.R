################################
### Kullback-Leibler divergence
################################

## band_choice = "DPI"

# kernel = "gaussian"

KLdiv_skew_t_DPI = matrix(NA, n_test, 2)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density_DPI[,ik], skew_t_DJI_ret_fore_den[,ik])
    colnames(skew_t_compar) = c("True", "Estimate")
    
    KLdiv_skew_t_DPI[ik,] = as.numeric(try(KLdiv(skew_t_compar, eps = 1e-16)), silent = TRUE)[2:3]
    rm(skew_t_compar)
}

KLdiv_skew_t_summary_DPI = round(sum(colMeans(KLdiv_skew_t_DPI, na.rm = TRUE)), 4) # 1.5841

# kernel = "epan"

KLdiv_skew_t_epan_DPI = matrix(NA, n_test, 2)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density_epan_DPI[,ik], skew_t_DJI_ret_fore_den[,ik])
    colnames(skew_t_compar) = c("True", "Estimate")
    
    KLdiv_skew_t_epan_DPI[ik,] = as.numeric(try(KLdiv(skew_t_compar, eps = 1e-16)), silent = TRUE)[2:3]
    rm(skew_t_compar)
}

KLdiv_skew_t_summary_epan_DPI = round(sum(colMeans(KLdiv_skew_t_epan_DPI, na.rm = TRUE)), 4) # 1.6668

##############################
### Jensen-Shannon divergence
##############################

## band_choice = "DPI"

# kernel = "gaussian"

JSdiv_skew_t_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density_DPI[,ik], skew_t_DJI_ret_fore_den[,ik])
    M = rowMeans(skew_t_compar)
    P_M = cbind(test_density_DPI[,ik], M)
    E_M = cbind(skew_t_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    A_M = try(KLdiv(P_M), silent = TRUE)
    B_M = try(KLdiv(E_M), silent = TRUE)
    if(class(A_M) == "try-error"|class(B_M) == "try-error")
    {
        JSdiv_skew_t_DPI[ik] = 0
    }
    else
    {
        JSdiv_skew_t_DPI[ik] = as.numeric(0.5 * A_M + 0.5 * B_M)[3]
    }
    rm(M); rm(P_M); rm(E_M); rm(skew_t_compar)
}

JSdiv_skew_t_summary_DPI = round(sum(JSdiv_skew_t_DPI[-which(JSdiv_skew_t_DPI == 0)]), 4) # 5.6072

# kernel = "epan"

JSdiv_skew_t_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density_epan_DPI[,ik], skew_t_DJI_ret_fore_den[,ik])
    M = rowMeans(skew_t_compar)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(skew_t_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    A_M = try(KLdiv(P_M), silent = TRUE)
    B_M = try(KLdiv(E_M), silent = TRUE)
    if(class(A_M) == "try-error"|class(B_M) == "try-error")
    {
        JSdiv_skew_t_epan_DPI[ik] = 0
    }
    else
    {
        JSdiv_skew_t_epan_DPI[ik] = as.numeric(0.5 * A_M + 0.5 * B_M)[3]
    }
    rm(M); rm(P_M); rm(E_M); rm(skew_t_compar)
}

JSdiv_skew_t_summary_epan_DPI = round(sum(JSdiv_skew_t_epan_DPI[-which(JSdiv_skew_t_epan_DPI == 0)]), 4) # 5.2563

########################################
# Jensen-Shannon divergence (geometric)
########################################

## band_choice = "DPI"

# kernel = "gaussian"

JSdiv_geo_skew_t_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density_DPI[,ik], skew_t_DJI_ret_fore_den[,ik])
    M = apply(skew_t_compar, 1, geometric.mean)
    P_M = cbind(test_density_DPI[,ik], M)
    E_M = cbind(skew_t_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    A_M = try(KLdiv(P_M), silent = TRUE)
    B_M = try(KLdiv(E_M), silent = TRUE)
    if(class(A_M) == "try-error"|class(B_M) == "try-error")
    {
        JSdiv_geo_skew_t_DPI[ik] = 0
    }
    else
    {
        JSdiv_geo_skew_t_DPI[ik] = as.numeric(0.5 * A_M + 0.5 * B_M)[3]
    }
    rm(M); rm(P_M); rm(E_M); rm(skew_t_compar)
}

JSdiv_geo_skew_t_DPI_summary = round(sum(JSdiv_geo_skew_t_DPI[-which(JSdiv_geo_skew_t_DPI == 0)]), 4) #

# kernel = "epan"

JSdiv_geo_skew_t_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density_epan_DPI[,ik], skew_t_DJI_ret_fore_den[,ik])
    M = apply(skew_t_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan_DPI[,ik], M)
    E_M = cbind(skew_t_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    A_M = try(KLdiv(P_M), silent = TRUE)
    B_M = try(KLdiv(E_M), silent = TRUE)
    if(class(A_M) == "try-error"|class(B_M) == "try-error")
    {
        JSdiv_geo_skew_t_epan_DPI[ik] = 0
    }
    else
    {
        JSdiv_geo_skew_t_epan_DPI[ik] = as.numeric(0.5 * A_M + 0.5 * B_M)[3]
    }
    rm(M); rm(P_M); rm(E_M); rm(skew_t_compar)
}

JSdiv_geo_skew_t_epan_DPI_summary = round(sum(JSdiv_geo_skew_t_epan_DPI[-which(JSdiv_geo_skew_t_epan_DPI == 0)]), 4) #


############
# L1 - norm
############

# kernel = "gaussian"

L1_norm_skew_t_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_skew_t_DPI[ik] = sum(abs(skew_t_DJI_ret_fore_den[,ik] - test_density_DPI[,ik]))
}

L1_norm_skew_t_summary_DPI = round(mean(L1_norm_skew_t_DPI[-30]), 4) # 1369.452

# kernel = "epan" (Another kernel function)

L1_norm_skew_t_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_skew_t_epan_DPI[ik] = sum(abs(skew_t_DJI_ret_fore_den[,ik] - test_density_epan_DPI[,ik]))
}

L1_norm_skew_t_summary_epan_DPI = round(mean(L1_norm_skew_t_epan_DPI[-30]), 4) # 1367.93

############
# L2 - norm
############

## band_choice = "DPI"

# kernel = "gaussian"

L2_norm_skew_t_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_skew_t_DPI[ik] = sqrt(sum((test_density_DPI[,ik] - skew_t_DJI_ret_fore_den[,ik])^2))
}

L2_norm_skew_t_summary_DPI = round(mean(L2_norm_skew_t_DPI[-30]), 4) # 67.9682

# kernel = "epan"

L2_norm_skew_t_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_skew_t_epan_DPI[ik] = sqrt(sum((test_density_epan_DPI[,ik] - skew_t_DJI_ret_fore_den[,ik])^2))
}

L2_norm_skew_t_summary_epan_DPI = round(mean(L2_norm_skew_t_epan_DPI[-30]), 4) # 63.7964

##############
# Linf - norm
##############

## band_choice = "DPI"

# kernel = "gaussian"

Linf_norm_skew_t_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_skew_t_DPI[ik] = max(abs(skew_t_DJI_ret_fore_den[,ik] - test_density_DPI[,ik]))
}

Linf_norm_skew_t_summary_DPI = round(mean(Linf_norm_skew_t_DPI[-30]), 4) # 5.7268

# kernel = "epan"

Linf_norm_skew_t_epan_DPI = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_skew_t_epan_DPI[ik] = max(abs(skew_t_DJI_ret_fore_den[,ik] - test_density_epan_DPI[,ik]))
}

Linf_norm_skew_t_summary_epan_DPI = round(mean(Linf_norm_skew_t_epan_DPI[-30]), 4) # 5.1539

##########
# summary
##########

DJI_skew_t_summary_DPI = c(KLdiv_skew_t_summary_DPI,
                           JSdiv_skew_t_summary_DPI,
                           JSdiv_geo_skew_t_DPI_summary,
                           L1_norm_skew_t_summary_DPI,
                           L2_norm_skew_t_summary_DPI,
                           Linf_norm_skew_t_summary_DPI)

DJI_skew_t_summary_epan_DPI = c(KLdiv_skew_t_summary_epan_DPI,
                                JSdiv_skew_t_summary_epan_DPI,
                                JSdiv_geo_skew_t_epan_DPI_summary,
                                L1_norm_skew_t_summary_epan_DPI,
                                L2_norm_skew_t_summary_epan_DPI,
                                Linf_norm_skew_t_summary_epan_DPI)

