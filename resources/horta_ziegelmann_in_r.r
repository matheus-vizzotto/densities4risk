library("ftsa")

# DADOS ORIGINAIS

data(DJI_return)

# Dow Jones Industrial Average (DJIA) is a stock market index that shows how 30 large publicly
# owned companies based in the United States have traded during a standard NYSE trading session.
# We consider monthly cross-sectional returns from April 2004 to December 2017. The data were
# obtained from the CRSP (Center for Research in Security Prices) database.

data <- DJI_return

nrow(data)
ncol(data)
DJI_return

Horta_Ziegelmann_FPCA(data = DJI_return, kernel = "epanechnikov",
    band_choice = "Silverman", ncomp_select = "FALSE")
# data: Densities or raw data matrix of dimension N by p, where N denotes sample size and p denotes dimensionality
# data: day x company -> Y: support x day

h_scale = 1
band_choice = "Silverman"
kernel = "epanechnikov"
m = 5001

    if(all(trunc(diff(apply(data, 1, sum))) == 0))
    {
        N = nrow(data)
        Y = t(data)
        u = gridpoints
        du = u[2] - u[1]
    } else  {
        # Sample size
        n = N = nrow(data)

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

        # defines gridpoints

        u = seq(from = min(data), to = max(data), length = m)
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
    }


# DADOS DE PESQUISA
library("readxl")
Y0 <- read_excel("C:/Users/user/Projetos/densities4risk/data/processed/lqdensities.xlsx")
Y <- as.matrix(Y0[ , -1])

Horta_Ziegelmann_FPCA(data = as.matrix(Y0[ , -1]), kernel = "epanechnikov",
    band_choice = "Silverman", ncomp_select = "FALSE")