# load libraries
library(RevGadgets)
library(psych)
library(coda)
library(dplyr)
library(ggplot2)

workingdir <- "/path/to/working/directory"
# Set working directory
setwd(workingdir)

### build index of sims where the harmonic mean of the ESS for parameters is > 200 
rel_sims <- c() # an empty vector for building the index of relevant sims
for (i in 1:100) {
  #load trace
    path <- paste0("sim_reduced/output/chisse_chisse/chisse_chisse_",i,".log") 
    trace <- readTrace(paths = path, burnin = 2000)
    trace <- trace[[1]]
    
    hmean <- harmonic.mean(coda::effectiveSize(trace)[2:ncol(trace)])
  
    if (hmean >= 200) { # if the harmonic ESS is >200, then add the ith sim to the growing index
      rel_sims <- rbind(rel_sims, i)
    }
} 

### calculate coverage, percent error in the mean, coefficient of variation in the posterior
# build data frame with mean, median, MAP, and 95% credible intervals for all parameters for each simulation
params <- c( "alpha", "chi",
             "clado_no_change_a", "clado_fission_a", "clado_fusion_a", "clado_demipoly_a", "clado_polyploid_a",
             "gamma_a", "delta_a", "eta_a", "rho_a",
             "clado_no_change_b", "clado_fission_b", "clado_fusion_b", "clado_demipoly_b", "clado_polyploid_b",
             "gamma_b", "delta_b", "eta_b", "rho_b",
             "turnover")
params_new <- c( "hidden_change_A", "hidden_change_C", 
             "no_change_C_i", "fission_C_i", "fusion_C_i", "demipoly_C_i", "polyploid_C_i",
             "fission_A_i", "fusion_A_i", "demipoly_A_i", "polyploid_A_i", 
             "no_change_C_ii", "fission_C_ii", "fusion_C_ii", "demipoly_C_ii", "polyploid_C_ii", 
             "fission_A_ii", "fusion_A_ii", "demipoly_A_ii", "polyploid_A_ii",
             "extinction")

df <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(df) <- c("sim_num", "parameter", "mean", "median", "MAP", "quantile_2.5", "quantile_97.5", "posterior_StDev")
for (i in rel_sims) { # for each simulation
  
  path <- paste0("sim_reduced/output/chisse_chisse/chisse_chisse_",i,".log") 
  trace <- readTrace(paths = path, burnin = 2000)
  # trace <- trace[[1]]
 
   for (j in 1:length(params)) { # within simulations for each parameter (so only have to load traces 1x)
     
     temp <- c()
     
     for (k in 1:5) { # within sims and parameters for each summary stat
       
       temp[k] <- summarizeTrace(trace, vars = params[j])[[1]][[1]][[k]]
       # temp <- c(1:5)
     }
    
     # st_dev <- sd(trace[,params[j]])
     st_dev <- sd(trace[[1]][,params[j]])
     
     df[(nrow(df)+1),] <- c(i, params[j], temp, st_dev)
  }
}
# df
# write.csv(df, "estimationError/parameterSummaries_raw.csv")

## Coverage
truth <- c(0.0103, 0.0242,
           0.582,    0.000671,0.00295, 0.00582, 0.0309,
           0.000671, 0.00537, 0.00492, 0.00277,
           0.0242,   0.0327,  0.0251,  0.0582,  0.224,
           0.336,    2.86,    0.121,   0.219,
           0.492) # these are the values underlying the simulation.

coverages <- data.frame(row.names = params)
coverages$perc_cover <- NA

for (i in 1:length(params)) {
  temp <- df[which(df$parameter==params[i]),]
  cover <- mean(truth[i] >= temp$quantile_2.5 & truth[i] <= temp$quantile_97.5)
  coverages$perc_cover[i] <- cover
}

# temp <- df %>%
#   group_by(parameter) %>%
#   summarise( mean_mean = mean(as.numeric(mean)) )

## % Error in the posterior MAP
percentError <- data.frame(row.names = params)
percentError$perc_err_mean <- NA
percentError$perc_err_median <- NA
percentError$perc_err_MAP <- NA
for (i in 1:length(params)) {
  temp <- df[which(df$parameter==params[i]),]
  perc_err_MAP <- ((as.numeric(temp$MAP) - truth[i]) / truth[i]) * 100
  perc_err_mean <- ((as.numeric(temp$mean) - truth[i]) / truth[i]) * 100
  perc_err_median <- ((as.numeric(temp$median) - truth[i]) / truth[i]) * 100
  percentError$perc_err_MAP[i] <- mean(perc_err_MAP)
  percentError$perc_err_mean[i] <- mean(perc_err_mean)
  percentError$perc_err_median[i] <- mean(perc_err_median)
}

# sapply(percentError, MARGIN = 2, "mean")

############# Plot % Error of the Mean for all params #############
library(ggplot2)
all_stErr <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(all_stErr) <- c("parameter","stErr")
for (i in 1:length(params)) {
  temp <- df[which(df$parameter==params[i]),]
  stErr <- ((as.numeric(temp$mean) - truth[i]) / truth[i]) * 100
  parameter <- rep(params[i], nrow(temp))
  all_stErr <- rbind(all_stErr, data.frame(parameter, stErr))
}
all_stErr$parameter <- factor(all_stErr$parameter, levels = params)
ggplot(all_stErr, aes(x = parameter, y = stErr)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_discrete(labels = params_new) +
  theme_classic() +
  labs( x = "Parameter", y = "Percent Error of the Mean (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
# boxplot(stErr ~ parameter, data = all_stErr)
###################################################################

## Coefficient of Variation of the posterior
CV <- data.frame(row.names = params)
CV$coeff_var <- NA

for (i in 1:length(params)) {
  temp <- df[which(df$parameter==params[i]),]
  coeff_var <- as.numeric(temp$posterior_StDev) / as.numeric(temp$mean)
  CV$coeff_var[i] <- mean(coeff_var)
}


## combine all these summaries
combined <- data.frame(row.names = params)
combined[,1:5] <- c(coverages$perc_cover, 
                    percentError$perc_err_mean, 
                    percentError$perc_err_median, 
                    percentError$perc_err_MAP, 
                    CV$coeff_var) %>% 
        round(digits = 2)

colnames(combined) <- c("perc_cover",
                        "perc_err_mean",
                        "perc_err_median",
                        "perc_err_MAP",
                        "coeff_var")


write.csv(combined, "estimationError/errorSummary.csv")


