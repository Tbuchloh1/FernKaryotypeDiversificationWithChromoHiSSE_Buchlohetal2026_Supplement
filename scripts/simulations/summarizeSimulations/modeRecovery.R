# load libraries
library(RevGadgets)
library(psych)
library(coda)
library(dplyr)
library(ggplot2)
library(ggbreak)

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



############# Plot estimates (absolute) for each parameter across simulations #############
library(ggplot2)
library(ggbreak)
ggplot(df, aes(x = as.factor(parameter), y = as.numeric(mean))) +
  geom_boxplot(outliers = F) +
  theme_classic() +
  labs( x = "Parameter", y = "Median Rate Estimate") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_break(c(0.12,0.3)) + 
  scale_y_break(c(0.5,0.8)) +
  ylim(0,1.2) +
  geom_vline(xintercept = seq(2.5, length(unique(df$parameter)) - 0.5, by = 2), 
             linetype = "dotted", color = "gray")


# boxplot(stErr ~ parameter, data = all_stErr)
