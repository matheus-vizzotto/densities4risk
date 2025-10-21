##########
# summary
##########

KLdiv_all_summary_gaussian = c(KLdiv_Horta_Ziegelman_summary, KLdiv_CoDa_summary, KLdiv_Alex_arima_summary, KLdiv_skew_t_summary)
JSdiv_all_summary_gaussian = c(JSdiv_Horta_Ziegelman_summary, JSdiv_CoDa_summary, JSdiv_Alex_arima_summary, JSdiv_skew_t_summary)
L1_norm_all_summary_gaussian = c(L1_norm_Horta_Ziegelman_summary, L1_norm_CoDa_summary, L1_norm_Alex_arima_summary, L1_norm_skew_t_summary)
L2_norm_all_summary_gaussian = c(L2_norm_Horta_Ziegelman_summary, L2_norm_CoDa_summary, L2_norm_Alex_arima_summary, L2_norm_skew_t_summary)
Linf_norm_all_summary_gaussian = c(Linf_norm_Horta_Ziegelman_summary, Linf_norm_CoDa_summary, Linf_norm_Alex_arima_summary, Linf_norm_skew_t_summary)
criteria_norm_all_summary_gaussian = rbind(KLdiv_all_summary_gaussian, JSdiv_all_summary_gaussian,
                                           L1_norm_all_summary_gaussian, L2_norm_all_summary_gaussian,
                                           Linf_norm_all_summary_gaussian)


KLdiv_all_summary_epan = c(KLdiv_Horta_Ziegelman_summary_epan, KLdiv_CoDa_summary_epan, KLdiv_Alex_arima_summary_epan, KLdiv_skew_t_summary_epan)
JSdiv_all_summary_epan = c(JSdiv_Horta_Ziegelman_summary_epan, JSdiv_CoDa_summary_epan, JSdiv_Alex_arima_summary_epan, JSdiv_skew_t_summary_epan)
L1_norm_all_summary_epan = c(L1_norm_Horta_Ziegelman_summary_epan, L1_norm_CoDa_summary_epan, L1_norm_Alex_arima_summary_epan, L1_norm_skew_t_summary_epan)
L2_norm_all_summary_epan = c(L2_norm_Horta_Ziegelman_summary_epan, L2_norm_CoDa_summary_epan, L2_norm_Alex_arima_summary_epan, L2_norm_skew_t_summary_epan)
Linf_norm_all_summary_epan = c(Linf_norm_Horta_Ziegelman_summary_epan, Linf_norm_CoDa_summary_epan, Linf_norm_Alex_arima_summary_epan, Linf_norm_skew_t_summary_epan)
criteria_norm_all_summary_epan = rbind(KLdiv_all_summary_epan, JSdiv_all_summary_epan,
                                      L1_norm_all_summary_epan, L2_norm_all_summary_epan,
                                      Linf_norm_all_summary_epan)
  
colnames(criteria_norm_all_summary_gaussian) = colnames(criteria_norm_all_summary_epan) = c("DFPCA", "CoDa", "LQD", "skew-t")
rownames(criteria_norm_all_summary_gaussian) = rownames(criteria_norm_all_summary_epan) = c("KL", "JS", "L1_norm", "L2_norm", "Linf_norm")

install.packages("xtable")
require(xtable)
xtable(criteria_norm_all_summary_gaussian, digits = 2)
xtable(criteria_norm_all_summary_epan, digits = 2)

######
# MCS
######

# Kullback-Leibler divergence

DJI_err_KLdiv = cbind(rowMeans(KLdiv_Horta_Ziegelman),
                      rowMeans(KLdiv_Alex_arima),
                      rowMeans(KLdiv_CoDa),
                      rowMeans(KLdiv_CoDa_no_normalization),
                      rowMeans(KLdiv_skew_t))

DJI_err_KLdiv_epan = cbind(rowMeans(KLdiv_Horta_Ziegelman_epan),
                           rowMeans(KLdiv_Alex_arima_epan),
                           rowMeans(KLdiv_CoDa_epan),
                           rowMeans(KLdiv_CoDa_epan_no_normalization),
                           rowMeans(KLdiv_skew_t_epan))

# Jensen-Shannon divergence (simple mean)

DJI_err_JSdiv = cbind(JSdiv_Horta_Ziegelman,
                      JSdiv_Alex_arima,
                      JSdiv_CoDa,
                      JSdiv_CoDa_no_normalization,
                      JSdiv_skew_t)

DJI_err_JSdiv_epan = cbind(JSdiv_Horta_Ziegelman_epan,
                           JSdiv_Alex_arima_epan,
                           JSdiv_CoDa_epan,
                           JSdiv_CoDa_epan_no_normalization,
                           JSdiv_skew_t_epan)

# Jensen-Shannon divergence (geometric mean)

DJI_err_JSdiv_geo = cbind(JSdiv_geo_Horta_Ziegelman,
                          JSdiv_geo_Alex_arima,
                          JSdiv_geo_CoDa,
                          JSdiv_geo_CoDa_no_normalization,
                          JSdiv_geo_skew_t)

DJI_err_JSdiv_geo_epan = cbind(JSdiv_geo_Horta_Ziegelman_epan,
                               JSdiv_geo_Alex_arima_epan,
                               JSdiv_geo_CoDa_epan,
                               JSdiv_geo_CoDa_epan_no_normalization,
                               JSdiv_geo_skew_t_epan)

# L1-norm

DJI_err_L1_norm = cbind(L1_norm_Horta_Ziegelman,
                        L1_norm_Alex_arima,
                        L1_norm_CoDa,
                        L1_norm_CoDa_no_normalization,
                        L1_norm_skew_t)

DJI_err_L1_norm_epan = cbind(L1_norm_Horta_Ziegelman_epan,
                             L1_norm_Alex_arima_epan,
                             L1_norm_CoDa_epan,
                             L1_norm_CoDa_epan_no_normalization,
                             L1_norm_skew_t_epan)

# L2-norm

DJI_err_L2_norm = cbind(L2_norm_Horta_Ziegelman,
                        L2_norm_Alex_arima,
                        L2_norm_CoDa,
                        L2_norm_CoDa_no_normalization,
                        L2_norm_skew_t)

DJI_err_L2_norm_epan = cbind(L2_norm_Horta_Ziegelman_epan,
                             L2_norm_Alex_arima_epan,
                             L2_norm_CoDa_epan,
                             L2_norm_CoDa_epan_no_normalization,
                             L2_norm_skew_t_epan)

# Linf-norm

DJI_err_Linf_norm = cbind(Linf_norm_Horta_Ziegelman,
                          Linf_norm_Alex_arima,
                          Linf_norm_CoDa,
                          Linf_norm_CoDa_no_normalization,
                          Linf_norm_skew_t)

DJI_err_Linf_norm_epan = cbind(Linf_norm_Horta_Ziegelman_epan,
                               Linf_norm_Alex_arima_epan,
                               Linf_norm_CoDa_epan,
                               Linf_norm_CoDa_epan_no_normalization,
                               Linf_norm_skew_t_epan)


# Kullback-Leibler divergence

DJI_MCS_95_KLdiv <- MCSprocedure(Loss = DJI_err_KLdiv[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 4
DJI_MCS_95_KLdiv_epan <- MCSprocedure(Loss = DJI_err_KLdiv_epan[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 4

# Jensen-Shannon divergence (simple mean)

DJI_MCS_95_JSdiv <- MCSprocedure(Loss = DJI_err_JSdiv[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 2
DJI_MCS_95_JSdiv_epan <- MCSprocedure(Loss = DJI_err_JSdiv_epan[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 2

# Jensen-Shannon divergence (geometric mean)

DJI_MCS_95_JSdiv_geo <- MCSprocedure(Loss = DJI_err_JSdiv_geo[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 4
DJI_MCS_95_JSdiv_geo_epan <- MCSprocedure(Loss = DJI_err_JSdiv_geo_epan[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 4

# L1-norm

DJI_MCS_95_L1_norm <- MCSprocedure(Loss = DJI_err_L1_norm[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 2 4
DJI_MCS_95_L1_norm_epan <- MCSprocedure(Loss = DJI_err_L1_norm_epan[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 2

# L2-norm

DJI_MCS_95_L2_norm <- MCSprocedure(Loss = DJI_err_L2_norm[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 2
DJI_MCS_95_L2_norm_epan <- MCSprocedure(Loss = DJI_err_L2_norm_epan[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 2

# Linf-norm

DJI_MCS_95_Linf_norm <- MCSprocedure(Loss = DJI_err_Linf_norm[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 4
DJI_MCS_95_Linf_norm_epan <- MCSprocedure(Loss = DJI_err_Linf_norm_epan[-30,], alpha = 0.95, B = 5000, statistic = "Tmax", cl = NULL) # 2

