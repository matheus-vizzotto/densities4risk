SmoLogFertRate = read.csv('SmoLogFertRate.csv')
n = dim(SmoLogFertRate)[2]

require(ftsa)

result_pivotal = T_stationary(sample = SmoLogFertRate, J = 100, 
                              MC_rep = 1000, pivotal = TRUE, Ker1=TRUE,
                              Ker2=FALSE)


