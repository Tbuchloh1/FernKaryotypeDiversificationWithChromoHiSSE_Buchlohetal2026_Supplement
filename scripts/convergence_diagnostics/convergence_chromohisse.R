library(coda)
library(RevGadgets)
library(psych)

setwd("/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/3.Results/chromohisse_v3")

# # POLYP
tracepaths <- c("output/polyp/polyp_chromohisse_tb_1.log",
                "output/polyp/polyp_chromohisse_tb_2.log",
                "output/polyp/polyp_chromohisse_tb_3.log",
                "output/polyp/polyp_chromohisse_tb_4.log",
                "output/polyp/polyp_chromohisse_tb_5.log")
# # ASPLE

# tracepaths <- c("output/asple/asple_chromohisse_tb_1.log",
#                 "output/asple/asple_chromohisse_tb_2.log",
#                 "output/asple/asple_chromohisse_tb_3.log",
#                 "output/asple/asple_chromohisse_tb_4.log",
#                 "output/asple/asple_chromohisse_tb_5.log")
# PTERI
# tracepaths <- c("output/pteri/pteri_chromohisse_tb_1.log",
#                 "output/pteri/pteri_chromohisse_tb_2.log",
#                 "output/pteri/pteri_chromohisse_tb_3.log",
#                 "output/pteri/pteri_chromohisse_tb_4.log",
#                 "output/pteri/pteri_chromohisse_tb_5.log")


# load mcmc log files
# test for convergence of the chains using ess and the chain's harmonic mean 
traces <- readTrace(tracepaths, burnin = 0.1)
traces[[6]] <- combineTraces(traces)[[1]]


hmeans <- numeric(length = length(traces))
for (i in 1:length(traces)) {
  hmeans[i] <- harmonic.mean(c(coda::effectiveSize(traces[[i]]))[2:ncol(traces[[i]])])
}
hmeans

#######
# TB
traces.df1 <- data.frame(readTrace(tracepaths[1], burnin = 0.1))
traces.df2 <- data.frame(readTrace(tracepaths[2], burnin = 0.1))
traces.df3 <- data.frame(readTrace(tracepaths[3], burnin = 0.1))
traces.df4 <- data.frame(readTrace(tracepaths[4], burnin = 0.1))
traces.df5 <- data.frame(readTrace(tracepaths[5], burnin = 0.1))

start = 0
end   = 12460
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
gelman.diag(mcmc_list, autoburnin = F)
gelman.plot(mcmc_list)
effectiveSize(mcmc_list)
lapply(mcmc_list, effectiveSize)

