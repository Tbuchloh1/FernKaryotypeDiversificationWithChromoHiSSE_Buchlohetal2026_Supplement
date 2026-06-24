library(coda)
library(RevGadgets)
library(psych)

setwd("/path/to/working/directory")

# # POLYP
# tracepaths <- c("output/polyp/polyp_chromosse_tb_1.log",
#                 "output/polyp/polyp_chromosse_tb_2.log",
#                 "output/polyp/polyp_chromosse_tb_3.log",
#                 "output/polyp/polyp_chromosse_tb_4.log",
#                 "output/polyp/polyp_chromosse_tb_5.log")
# ASPLE
# tracepaths <- c("output/asple/asple_chromosse_tb_1.log",
#                 "output/asple/asple_chromosse_tb_2.log",
#                 "output/asple/asple_chromosse_tb_3.log",
#                 "output/asple/asple_chromosse_tb_4.log",
#                 "output/asple/asple_chromosse_tb_5.log")
# PTERI
tracepaths <- c("output/pteri/pteri_chromosse_tb_1.log",
                "output/pteri/pteri_chromosse_tb_2.log",
                "output/pteri/pteri_chromosse_tb_3.log",
                "output/pteri/pteri_chromosse_tb_4.log",
                "output/pteri/pteri_chromosse_tb_5.log")

# # Sims
# setwd("/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/2.Methods/tests/sim_reduced/r1_pteri/output/csse_csse")
# tracepaths <- "csse_csse_101.log"


# # load mcmc log files
# # test for convergence of the chains using ess and the chain's harmonic mean 
# traces <- readTrace(tracepaths, burnin = 2000)
# traces[[6]] <- combineTraces(traces)[[1]]
# 
# 
# hmeans <- numeric(length = length(traces))
# for (i in 1:length(traces)) {
#   hmeans[i] <- harmonic.mean(c(coda::effectiveSize(traces[[i]]))[2:ncol(traces[[i]])])
# }
# hmeans

#######
# diagnostics 
vars <- c("Likelihood",
          "clado_demipoly", "clado_fission", "clado_fusion", "clado_no_change", "clado_polyploid",
          "delta", "eta", "gamma", "rho", 
          "turnover")
traces.df1 <- data.frame(readTrace(tracepaths[1], burnin = 0.1))
traces.df1 <- traces.df1[, vars]
traces.df2 <- data.frame(readTrace(tracepaths[2], burnin = 0.1))
traces.df2 <- traces.df2[, vars]
traces.df3 <- data.frame(readTrace(tracepaths[3], burnin = 0.1))
traces.df3 <- traces.df3[, vars]
traces.df4 <- data.frame(readTrace(tracepaths[4], burnin = 0.1))
traces.df4 <- traces.df4[, vars]
traces.df5 <- data.frame(readTrace(tracepaths[5], burnin = 0.1))
traces.df5 <- traces.df5[, vars]

start = 0
end   = min(c(nrow(traces.df1),nrow(traces.df2),nrow(traces.df3),nrow(traces.df4),nrow(traces.df5))) - 1
traces.mcmc1 <- mcmc(traces.df1, start = start, end = end)
traces.mcmc2 <- mcmc(traces.df2, start = start, end = end)
traces.mcmc3 <- mcmc(traces.df3, start = start, end = end)
traces.mcmc4 <- mcmc(traces.df4, start = start, end = end)
traces.mcmc5 <- mcmc(traces.df5, start = start, end = end)

mcmc_list <- mcmc.list(traces.mcmc1,
                       traces.mcmc2,
                       traces.mcmc3,
                       traces.mcmc4,
                       traces.mcmc5)

# gelman.diag(mcmc_list, autoburnin = F)
# gelman.plot(mcmc_list)
# effectiveSize(mcmc_list)
# lapply(mcmc_list, effectiveSize)

##########
####### 
## catch outputs here one at a time

# make dataframe to catch diagnostics 
# chromosse_mcmcdiagnostics <- data.frame(vars = vars, polyp_ess = NA, polyp_rhat_pt = NA, polyp_rhat_uCI = NA,
# asple_ess = NA, asple_rhat_pt = NA, asple_rhat_uCI = NA,
# pteri_ess = NA, pteri_rhat_pt = NA, pteri_rhat_uCI = NA)

# polypodiineae
# chromosse_mcmcdiagnostics$polyp_ess <- effectiveSize(mcmc_list)
# chromosse_mcmcdiagnostics$polyp_rhat_pt <- gelman.diag(mcmc_list, autoburnin = F)$psrf[,1]
# chromosse_mcmcdiagnostics$polyp_rhat_uCI <- gelman.diag(mcmc_list, autoburnin = F)$psrf[,2]

# aspleniineae
# chromosse_mcmcdiagnostics$asple_ess <- effectiveSize(mcmc_list)
# chromosse_mcmcdiagnostics$asple_rhat_pt <- gelman.diag(mcmc_list, autoburnin = F)$psrf[,1]
# chromosse_mcmcdiagnostics$asple_rhat_uCI <- gelman.diag(mcmc_list, autoburnin = F)$psrf[,2]

# pteridineae
chromosse_mcmcdiagnostics$pteri_ess <- effectiveSize(mcmc_list)
chromosse_mcmcdiagnostics$pteri_rhat_pt <- gelman.diag(mcmc_list, autoburnin = F)$psrf[,1]
chromosse_mcmcdiagnostics$pteri_rhat_uCI <- gelman.diag(mcmc_list, autoburnin = F)$psrf[,2]

#######
# write diagnostics table to .csv
write.csv(chromosse_mcmcdiagnostics, "/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/2.Methods/postprocessing/convergence/chromosse_mcmcconvergence_diagnostics.csv")
