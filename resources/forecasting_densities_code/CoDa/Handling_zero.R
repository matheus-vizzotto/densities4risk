###########################
# Dealing with zero values
###########################

# c = 10^6

DJI_return_density_transformation = DJI_return_density_trans * (10^6)

# maximum rounding-off error

epsilon = sapply(1:165, function(X) max(DJI_return_density_transformation[,X] - round(DJI_return_density_transformation[,X], 2)))

# multiplicative replacement strategy

DJI_CoDa_mat = matrix(NA, 1000, 165)
for(ik in 1:165)
{ 
    index = which(round(DJI_return_density_transformation[,ik], 2) == 0)
    DJI_CoDa_mat[,ik] = replace(DJI_return_density_transformation[,ik], index, epsilon[ik])
    DJI_CoDa_mat[-index,ik] = DJI_return_density_transformation[-index,ik] * (1 - (length(index) * epsilon[ik])/(10^6))
    print(ik)
}    
