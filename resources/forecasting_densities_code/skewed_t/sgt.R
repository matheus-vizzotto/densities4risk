#################
# load R package
#################

install.packages("fGarch")
require(fGarch)
require(vars)

# Importing the data
setwd("/Volumes/MY PASSPORT/Dropbox/Todos/cs/data")
DJI_data = read.table("DJI_monthly.txt", header = TRUE, sep = "")
DJI_ret = as.matrix(DJI_data[,2:31])


skew_t_fun <- function(data, M = 5001)
{
    n = nrow(data)
    
    # Parameter estimation

    para_est = sapply(1:n, function(ik) sstdFit(data[ik,])$estimate)
    m = seq(min(data), max(data), length.out = M)
    den_skewed_t = matrix(NA, M, n)
    for(iw in 1:n)
    {
        den_skewed_t[,iw] = dsstd(x = m, mean = para_est[1,iw], sd = para_est[2,iw], 
                                  nu = para_est[3,iw], xi = para_est[4,iw])
        print(iw)
    }

    # Re-arrange parameters

    para_est_new = matrix(NA, 4, n)

    # log of variances

    para_est_new[2,] = log(para_est[2,]^2)

    # log of transformed df

    x = vector("numeric", length(para_est[3,]))
    for(i in 1:length(para_est[3,]))
    {
        x[i] = pt(q = -2, df = para_est[3,i])
    }
    para_est_new[3,] = log(x)

    # mean and skewness
    
    para_est_new[1,] = para_est[1,]
    para_est_new[4,] = para_est[4,]

    # differencing log of variance and differencing log of transformed df
    
    parest_skewt_chen = matrix(nr = 4, nc = ncol(para_est) - 1)
    parest_skewt_chen[2,] = diff(para_est_new[2,])
    parest_skewt_chen[3,] = diff(para_est_new[3,])
    parest_skewt_chen[1,] = para_est_new[1,-1]
    parest_skewt_chen[4,] = para_est_new[4,-1]
    rownames(parest_skewt_chen) = c("mean", "diff(log(variance))", "diff(log(transformation))", "skewness")

    # Back-transformation of log of transformed df
    
    A = exp(para_est_new[3,])
    df_cal = vector()
    for(i in 1:length(A))
    {
        df_cal[i] = uniroot(f = function(df){pt(-2, df) - A[i]}, 
                            interval = c(min(para_est[3,])-1, max(para_est[3,])+1))$root
    }

    # forecasting four parameters

    VAR_order = VARselect(t(parest_skewt_chen))$selection[1]
    VAR_fore = predict(VAR(t(parest_skewt_chen), p = VAR_order), n.ahead = 1)
    VAR_fore_parest = sapply(1:4, function(ik) VAR_fore$fcst[[ik]][1])
    
    VAR_fore_parest_new = vector("numeric", 4)
    VAR_fore_parest_new[2] = sqrt(exp(VAR_fore_parest[2] + tail(para_est_new[2,], 1)))
    
    A_fore = exp(VAR_fore_parest[3] + tail(para_est_new[3,],1))
    dum = try(uniroot(f = function(df){pt(-2, df) - A_fore}, 
                          interval = c(min(para_est[3,])-1, max(para_est[3,])+1))$root, silent = TRUE)
    if(class(dum) == "try-error")
    {
        print("Unit root problem")   
        VAR_fore_parest_new[3] = tail(df_cal,1)
    }
    else
    {
        VAR_fore_parest_new[3] = dum
    }
    VAR_fore_parest_new[1] = VAR_fore_parest[1]
    VAR_fore_parest_new[4] = VAR_fore_parest[4]

    skewed_t_den_fore = dsstd(x = m, mean = VAR_fore_parest_new[1], sd = VAR_fore_parest_new[2], 
                              nu = VAR_fore_parest_new[3], xi = VAR_fore_parest_new[4])
    return(list(m = m, skewed_t_den_fore = skewed_t_den_fore))
}

# forecasting 

skew_t_DJI_ret_fore_den = matrix(NA, 5001, 55)
for(ik in 1:55)
{
    skew_t_DJI_ret_fore_den[,ik] = skew_t_fun(data = DJI_ret[1:(109+ik),])$skewed_t_den_fore
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
test_density_epan = ken_density(data = DJI_ret, band_choice = "Silverman", kernel = "epanechnikov")$Y[,111:165]
n_test = ncol(test_density_epan)

test_density_DPI = ken_density(data = DJI_ret, band_choice = "DPI", kernel = "gaussian")$Y[,111:165]
test_density_epan_DPI = ken_density(data = DJI_ret, band_choice = "DPI", kernel = "epan")$Y[,111:165]

#####################
# Density evaluation
#####################

### Kullback-Leibler divergence

## band_choice = "Silverman"

# kernel = "gaussian"

KLdiv_skew_t = matrix(NA, n_test, 2)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density[,ik], skew_t_DJI_ret_fore_den[,ik])
    colnames(skew_t_compar) = c("True", "Estimate")
  
    dum = try(as.numeric(KLdiv(skew_t_compar, eps = 1e-16))[2:3], silent = TRUE)
    if(class(dum) == "try-error")
    {
        KLdiv_skew_t[ik,] = c(NA, NA)
    }
    else
    {
        KLdiv_skew_t[ik,] = dum
    }
    print(ik); rm(skew_t_compar); rm(dum)
}

KLdiv_skew_t_summary = round(sum(colMeans(KLdiv_skew_t, na.rm = TRUE)), 4)  # 1.359

# kernel = "epan"

KLdiv_skew_t_epan  = matrix(NA, n_test, 2)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density_epan[,ik], skew_t_DJI_ret_fore_den[,ik])
    colnames(skew_t_compar) = c("True", "Estimate")
    
    dum = try(as.numeric(KLdiv(skew_t_compar, eps = 1e-16))[2:3], silent = TRUE)
    if(class(dum) == "try-error")
    {
        KLdiv_skew_t_epan[ik,] = c(NA, NA)
    }
    else
    {
        KLdiv_skew_t_epan[ik,] = dum
    }
    print(ik); rm(skew_t_compar); rm(dum)
}

KLdiv_skew_t_summary_epan = round(sum(colMeans(KLdiv_skew_t_epan, na.rm = TRUE)), 4) # 1.5595

##########################################
# Jensen-Shannon divergence (simple mean)
##########################################

## band_choice = "Silverman"

# kernel = "gaussian"

JSdiv_skew_t = vector("numeric", n_test)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density[,ik], skew_t_DJI_ret_fore_den[,ik])
    M = rowMeans(skew_t_compar)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(skew_t_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    A_M = try(KLdiv(P_M), silent = TRUE)
    B_M = try(KLdiv(E_M), silent = TRUE)
    if(class(A_M) == "try-error"|class(B_M) == "try-error")
    {
        JSdiv_skew_t[ik] = 0    
    }
    else
    {
        JSdiv_skew_t[ik] = as.numeric(0.5 * A_M + 0.5 * B_M)[3]
    }
    rm(M); rm(P_M); rm(E_M); rm(skew_t_compar)
}

JSdiv_skew_t_summary = round(sum(JSdiv_skew_t[-which(JSdiv_skew_t == 0)]), 4) # 5.2532 

# kernel = "epan"

JSdiv_skew_t_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density_epan[,ik], skew_t_DJI_ret_fore_den[,ik])
    M = rowMeans(skew_t_compar)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(skew_t_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    A_M = try(KLdiv(P_M), silent = TRUE)
    B_M = try(KLdiv(E_M), silent = TRUE)
    if(class(A_M) == "try-error"|class(B_M) == "try-error")
    {
        JSdiv_skew_t_epan[ik] = 0
    }
    else
    {
        JSdiv_skew_t_epan[ik] = as.numeric(0.5 * A_M + 0.5 * B_M)[3]
    }
    rm(M); rm(P_M); rm(E_M); rm(skew_t_compar)
}

JSdiv_skew_t_summary_epan = round(sum(JSdiv_skew_t_epan[-which(JSdiv_skew_t_epan == 0)]), 4) # 5.4505

########################################
# Jensen-Shannon divergence (geometric)
########################################

## band_choice = "Silverman"

# kernel = "gaussian"

JSdiv_geo_skew_t = vector("numeric", n_test)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density[,ik], skew_t_DJI_ret_fore_den[,ik])
    M = apply(skew_t_compar, 1, geometric.mean)
    P_M = cbind(test_density[,ik], M)
    E_M = cbind(skew_t_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    A_M = try(KLdiv(P_M), silent = TRUE)
    B_M = try(KLdiv(E_M), silent = TRUE)
    if(class(A_M) == "try-error"|class(B_M) == "try-error")
    {
        JSdiv_geo_skew_t[ik] = 0    
    }
    else
    {
        JSdiv_geo_skew_t[ik] = as.numeric(0.5 * A_M + 0.5 * B_M)[3]
    }
    rm(M); rm(P_M); rm(E_M); rm(skew_t_compar)
}

JSdiv_geo_skew_t_summary = round(sum(JSdiv_geo_skew_t[-which(JSdiv_geo_skew_t == 0)]), 4) # 5.2532 

# kernel = "epan"

JSdiv_geo_skew_t_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    skew_t_compar = cbind(test_density_epan[,ik], skew_t_DJI_ret_fore_den[,ik])
    M = apply(skew_t_compar, 1, geometric.mean)
    P_M = cbind(test_density_epan[,ik], M)
    E_M = cbind(skew_t_DJI_ret_fore_den[,ik], M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    A_M = try(KLdiv(P_M), silent = TRUE)
    B_M = try(KLdiv(E_M), silent = TRUE)
    if(class(A_M) == "try-error"|class(B_M) == "try-error")
    {
        JSdiv_geo_skew_t_epan[ik] = 0
    }
    else
    {
        JSdiv_geo_skew_t_epan[ik] = as.numeric(0.5 * A_M + 0.5 * B_M)[3]
    }
    rm(M); rm(P_M); rm(E_M); rm(skew_t_compar)
}

JSdiv_geo_skew_t_summary_epan = round(sum(JSdiv_geo_skew_t_epan[-which(JSdiv_geo_skew_t_epan == 0)]), 4) # 5.4505


############
# L1 - norm
############

## band_choice = "Silverman"

# kernel = "gaussian"

L1_norm_skew_t = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_skew_t[ik] = sum(abs(skew_t_DJI_ret_fore_den[,ik] - test_density[,ik]))
}  

L1_norm_skew_t_summary = round(mean(L1_norm_skew_t[-30]), 4) # 1324.9720

# kernel = "epan"

L1_norm_skew_t_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L1_norm_skew_t_epan[ik] = sum(abs(skew_t_DJI_ret_fore_den[,ik] - test_density_epan[,ik]))
}

L1_norm_skew_t_summary_epan = round(mean(L1_norm_skew_t_epan[-30]), 4) # 1423.941


##########
# L2-norm
##########

## band_choice = "Silverman"

# kernel = "gaussian"

L2_norm_skew_t = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_skew_t[ik] = sqrt(sum((test_density[,ik] - skew_t_DJI_ret_fore_den[,ik])^2))
}

L2_norm_skew_t_summary = round(mean(L2_norm_skew_t[-30]), 4) # 64.3983

# kernel = "epan"

L2_norm_skew_t_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    L2_norm_skew_t_epan[ik] = sqrt(sum((test_density_epan[,ik] - skew_t_DJI_ret_fore_den[,ik])^2))    
}

L2_norm_skew_t_summary_epan = round(mean(L2_norm_skew_t_epan[-30]), 4) # 65.6475


##############
# Linf - norm
##############

## band_choice = "Silverman"

# kernel = "gaussian"

Linf_norm_skew_t = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_skew_t[ik] = max(abs(skew_t_DJI_ret_fore_den[,ik] - test_density[,ik]))
}  

Linf_norm_skew_t_summary = round(mean(Linf_norm_skew_t[-30]), 4) # 5.1945

# kernel = "epan"

Linf_norm_skew_t_epan = vector("numeric", n_test)
for(ik in 1:n_test)
{
    Linf_norm_skew_t_epan[ik] = max(abs(skew_t_DJI_ret_fore_den[,ik] - test_density_epan[,ik]))
}

Linf_norm_skew_t_summary_epan = round(mean(Linf_norm_skew_t_epan[-30]), 4) # 5.38

##########
# summary
##########

DJI_skew_t_summary = c(KLdiv_skew_t_summary,
                       JSdiv_skew_t_summary,
                       JSdiv_geo_skew_t_summary,
                       L1_norm_skew_t_summary,
                       L2_norm_skew_t_summary,
                       Linf_norm_skew_t_summary)

DJI_skew_t_summary_epan = c(KLdiv_skew_t_summary_epan,
                            JSdiv_skew_t_summary_epan,
                            JSdiv_geo_skew_t_summary_epan,
                            L1_norm_skew_t_summary_epan,
                            L2_norm_skew_t_summary_epan,
                            Linf_norm_skew_t_summary_epan)

