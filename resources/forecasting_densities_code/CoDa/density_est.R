###############
# 30 companies
###############

setwd("/Volumes/MY PASSPORT/Dropbox/Todos/cs/data")
DJI = read.table("DJI_monthly.txt", header = TRUE)
DJI_month = DJI[,1]
DJI_return = as.matrix(DJI[,2:31])

#######################
# Exploratory analysis
#######################

kde_DJI <- function(data, m)
{
    N = nrow(data)
    h_scale = 1
    h.hat_5m = sapply(1:N, function(t) 1.06*sd(data[t,])*(length(data[t,])^(-(1/5))))

    h.hat_5m = h_scale * h.hat_5m

    # Initialization parameters
    n = N # Number of daily observations
    
    # 2. Discretization
    # Evaluation points
    u = seq(from = min(data), to = max(data), length = m)
    
    # Interval length
    du = u[2] - u[1]
    
    Y = sapply(1:N, function(t) density(data[t,], bw = h.hat_5m[t], kernel = 'gaussian', from = min(data), to = max(data), n = m)$y)
    # correcting to ensure integral Y_t du = 1
    for(t in 1:N)
    {
        Y[,t] = Y[,t]/(sum(Y[,t])*du)
    }
    return(list(Y = Y, u = u))
}    

expo_kde_DJI = kde_DJI(data = DJI_return, m = 5001)
persp(x = expo_kde_DJI$u, y = 1:165, z = expo_kde_DJI$Y, xlab = "Grid point", ylab = "Kernel density estimate")

# plot conditional hdrcde (play round)

plot(cde(x = expo_kde_DJI$u, y = expo_kde_DJI$Y[,1]), xlab = "Grid point", ylab = "Kernel density estimate")


month_index = c(which(DJI_month == "1/7/08"), which(DJI_month == "1/10/08"),
                which(DJI_month == "1/10/11"), which(DJI_month == "1/10/15"),
                which(DJI_month == "1/5/17"), which(DJI_month == "1/8/17"))

savepdf("DJI_fig_all", width = 12, height = 10, toplines = 0.8)
plot(expo_kde_DJI$u, expo_kde_DJI$Y[,month_index[2]], xlab = "Grid point", ylim = c(0, 11), 
     ylab = "Density", type = "l", lty = 1, col = 1)
lines(expo_kde_DJI$u, expo_kde_DJI$Y[,month_index[4]], lty = 2, col = 2)
lines(expo_kde_DJI$u, expo_kde_DJI$Y[,month_index[6]], lty = 4, col = 4)
legend("topleft", c("1/10/08", "1/10/15", "1/08/17"), col = c(1,2,4), lty = c(1,2,4), cex = 0.8)
dev.off()

#########################
# Plot graphs separately
#########################

# Negative returns

setwd("/Volumes/MY PASSPORT/Dropbox/Todos/cs/tex/plots")
savepdf("DJI_fig_1_a", width = 12, height = 10, toplines = 0.8)
plot(expo_kde_DJI$u, expo_kde_DJI$Y[,month_index[1]], xlab = "Grid point", ylab = "Density", type = "l", main = "1/July/2008")
dev.off()

savepdf("DJI_fig_1_b", width = 12, height = 10, toplines = 0.8)
plot(expo_kde_DJI$u, expo_kde_DJI$Y[,month_index[2]], xlab = "Grid point", ylab = "Density", type = "l", main = "1/Oct/2008")
dev.off()

# Positive returns

savepdf("DJI_fig_2_a", width = 12, height = 10, toplines = 0.8)
plot(expo_kde_DJI$u, expo_kde_DJI$Y[,month_index[3]], xlab = "Grid point", ylab = "Density", type = "l", main = "1/Oct/2011")
dev.off()

savepdf("DJI_fig_2_b", width = 12, height = 10, toplines = 0.8)
plot(expo_kde_DJI$u, expo_kde_DJI$Y[,month_index[4]], xlab = "Grid point", ylab = "Density", type = "l", main = "1/Oct/2015")
dev.off()

# Zero returns

savepdf("DJI_fig_3_a", width = 12, height = 10, toplines = 0.8)
plot(expo_kde_DJI$u, expo_kde_DJI$Y[,month_index[5]], xlab = "Grid point", ylab = "Density", type = "l", main = "1/May/2017")
dev.off()

savepdf("DJI_fig_3_b", width = 12, height = 10, toplines = 0.8)
plot(expo_kde_DJI$u, expo_kde_DJI$Y[,month_index[6]], xlab = "Grid point", ylab = "Density", type = "l", main = "1/Aug/2017")
dev.off()

