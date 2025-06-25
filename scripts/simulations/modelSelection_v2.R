library(RevGadgets)
library(geiger)
library(coda)
library(psych)


workingdir <- "/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/3.Results/sims/sim_reduced"
# Set working directory
setwd(workingdir)

# create empty dataframe
df <- data.frame(matrix(nrow = 200, ncol = 23))
colnames(df) <- c("sim", "gen_model", 
                  "hisse_lik", "hisse_lik_ess", 
                  "hisse_harmonic_ess", "sse_lik",  "sse_lik_ess",
                  "hisse_AIC", "hisse_dAIC", "sse_AIC", "sse_dAIC", 
                  "hisse_AICc", "hisse_dAICc", "sse_AICc", "sse_dAICc", 
                  "hisse_AICM", "hisse_dAICM","sse_AICM", "sse_dAICM",
                  "hisse_BIC", "hisse_dBIC", "sse_BIC", "sse_dBIC")

# label sims
df$sim <- 1:200
df$gen_model[1:100] <- "ChromoHiSSE"
df$gen_model[101:200] <- "ChromoSSE"

# loop through traces and fill out table. 
for (i in c(1:155, 157:200)) { # skip 156 since the analysis had weird errors. 
  #load trace
 if (i <= 100) {
    chisse_path <- paste0("output/chisse_chisse/chisse_chisse_",i,".log") 
    trace_1 <- readTrace(paths = chisse_path, burnin = 2000)
    trace_1 <- trace_1[[1]]
    
    csse_path <- paste0("output/chisse_csse/chisse_csse_",i,".log") 
    trace_2 <- readTrace(paths = csse_path, burnin = 2000)
    trace_2 <- trace_2[[1]]
    
    data_path <- paste0("data/chisse/csse/sim",i,".tsv") 
    
    df$hisse_harmonic_ess[i] <- harmonic.mean(coda::effectiveSize(trace_1)[2:ncol(trace_1)])
    
    } else if (i > 100) {
    chisse_path <- paste0("output/csse_chisse/csse_chisse_",i,".log") 
    trace_1 <- readTrace(paths = chisse_path, burnin = 2000)
    trace_1 <- trace_1[[1]]
    
    csse_path <- paste0("output/csse_csse/csse_csse_",i,".log") 
    trace_2 <- readTrace(paths = csse_path, burnin = 2000)
    trace_2 <- trace_2[[1]]
    
    data_path <- paste0("data/csse/csse/sim",i,".tsv")
    }
  
  if (i %in% c(15, 68, 87, 92)) { # catch weird analyses that got stuck and prune more of the burnin
    chisse_path <- paste0("output/chisse_chisse/chisse_chisse_",i,".log") 
    trace_1 <- readTrace(paths = chisse_path, burnin = 2000)
    trace_1 <- trace_1[[1]]
    
    csse_path <- paste0("output/chisse_csse/chisse_csse_",i,".log") 
    trace_2 <- readTrace(paths = csse_path, burnin = 15000)
    trace_2 <- trace_2[[1]]
    
    data_path <- paste0("data/chisse/csse/sim",i,".tsv") 
    
    df$hisse_harmonic_ess[i] <- harmonic.mean(coda::effectiveSize(trace_1)[2:ncol(trace_1)])
    
  } else if (i == 113) { # catch weird analyses that got stuck and prune more of the burnin
    chisse_path <- paste0("output/csse_chisse/csse_chisse_",i,".log") 
    trace_1 <- readTrace(paths = chisse_path, burnin = 2000)
    trace_1 <- trace_1[[1]]
    
    csse_path <- paste0("output/csse_csse/csse_csse_",i,".log") 
    trace_2 <- readTrace(paths = csse_path, burnin = 12200)
    trace_2 <- trace_2[[1]]
    
    data_path <- paste0("data/csse/csse/sim",i,".tsv")
  } 
  
  taxa <- read.table(data_path, header = F)
  n_taxa <- nrow(taxa)
  
  #calculate summary stats
  L1 <- getMAP(trace_1$Likelihood)
  L2 <- getMAP(trace_2$Likelihood)
  k1 <- 21
  k2 <- 10
  
  
  # fill in general info
  df$hisse_lik[i] <- L1
  df$hisse_lik_ess[i] <- effectiveSize(trace_1$Likelihood)
  
  df$sse_lik[i] <- L2
  df$sse_lik_ess[i] <- effectiveSize(trace_2$Likelihood)
  
  # calculate information criteria
  df$hisse_AIC[i] <- 2 * k1 - 2 * L1 
  df$sse_AIC[i] <- 2 * k2 - 2 * L2
  
  df$hisse_AICc[i] <- (2 * k1 - 2 * L1) + (2 * k1 * (k1 + 1)) / (n_taxa - k1 - 1)
  df$sse_AICc[i] <- (2 * k2 - 2 * L2) + (2 * k2 * (k2 + 1)) / (n_taxa - k2 - 1)
  
  df$hisse_AICM[i] <- aicm(trace_1$Likelihood[seq(1,nrow(trace_1), by = 50)])
  df$sse_AICM[i] <- aicm(trace_2$Likelihood[seq(1,nrow(trace_2), by = 50)])
  
  df$hisse_BIC[i] <- log(n_taxa) * k1 - 2 * L1
  df$sse_BIC[i] <- log(n_taxa) * k2 - 2 * L2
  
  if (i %% 10 == 0) {
    print(i)
  }

}

x = df$hisse_BIC[1:100]
y = df$hisse_BIC[101:200]
### model selection
# calculate dAICc
df[, c(9,11)] <- t(apply(df[, c(8,10)], 1, function(x) x - min(x)))
# calculate dAICc
df[, c(13,15)] <- t(apply(df[, c(12,14)], 1, function(x) x - min(x)))
# calculate dAICM
df[, c(17,19)] <- t(apply(df[, c(16,18)], 1, function(x) x - min(x)))
# calculate dBIC
df[, c(21,23)] <- t(apply(df[, c(20,22)], 1, function(x) x - min(x)))

### Summarize model selection
# get sims where likelihood ess is sufficient
hisse_sel <- df[which(df$hisse_lik_ess[1:100]>=200 & df$sse_lik_ess[1:100]>=200),c(9,13,17,21,11,15,19,23)] # 96 sims
sse_sel <- df[(which(df$hisse_lik_ess[101:200]>=200 & df$sse_lik_ess[101:200]>=200)+100),c(9,13,17,21,11,15,19,23)] # 98 sims

# make summary table
modelAd <- data.frame(matrix(nrow=2, ncol=9))
colnames(modelAd) <- c("Gen_model","dAIC = 0","dAICc = 0","dAICM = 0","dBIC = 0","dAIC > 10","dAICc > 10","dAICM > 10","dBIC > 10")
modelAd[,1] <- c("ChromoHiSSE", "ChromoSSE")
modelAd[1,2:5] <- apply(hisse_sel[,1:4], 2, function(row) {
  mean(row == 0, na.rm = TRUE)
})
modelAd[1,6:9] <- apply(hisse_sel[,5:8], 2, function(row) {
  mean(row >= 10)
})
modelAd[2,2:5] <- apply(sse_sel[,5:8], 2, function(row) {
  mean(row == 0, na.rm = TRUE)
})
modelAd[2,6:9] <- apply(sse_sel[,1:4], 2, function(row) {
  mean(row >= 10)
})
modelAd[,2:9] <- round(modelAd[,2:9], digits = 2)
modelAd

length(which(df$hisse_harmonic_ess>=200))

###### Write these data to tables and .csv so I can work with other params and make a table for the manuscript ########
# Write adequacy table
# write.csv(modelAd, "/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/3.Results/sims/modelAdequacy_sims.csv")

# Write dataframe with information criteria and ESS
# write.csv(df, "/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/3.Results/sims/rawModSel_sims.csv")

