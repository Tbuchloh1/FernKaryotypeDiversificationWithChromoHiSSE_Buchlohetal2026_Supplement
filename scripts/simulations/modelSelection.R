library(RevGadgets)

workingdir <- "/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/2.Methods/tests/sim_reduced"
# Set working directory
setwd(workingdir)

# create empty dataframe
df <- data.frame(matrix(nrow = 100, ncol = 15))
colnames(df) <- c("sim", 
                  "chisse_lnL", "csse_lnL", "chisse_AIC", "chisse_dAIC", "csse_AIC", "csse_dAIC", 
                  "chisse_AICc", "chisse_dAICc", "csse_AICc", "csse_dAICc", 
                  "chisse_BIC", "chisse_dBIC", "csse_BIC", "csse_dBIC")
# df$sim <- 1:100
df$sim <- 101:200

for (i in 1:100) {
  #load trace
  # chisse_path <- paste0("output/chisse_chisse/chisse_chisse_",i,".log")
  # csse_path <- paste0("output/chisse_csse/chisse_csse_",i,".log")
  chisse_path <- paste0("output/csse_chisse/csse_chisse_",i+100,".log")
  csse_path <- paste0("output/csse_csse/csse_csse_",i+100,".log")
  trace_1 <- readTrace(paths = chisse_path, burnin = 2000)
  trace_1 <- trace_1[[1]]
  trace_2 <- readTrace(paths = csse_path, burnin = 2000)
  trace_2 <- trace_2[[1]]
  
  #load dataset for calculating number of taxa in the dataset
  # data_path <- paste0("data/chisse/csse/sim",i,".tsv")
  data_path <- paste0("data/csse/csse/sim",i+100,".tsv")
  taxa <- read.table(data_path, header = F)
  n_taxa <- nrow(taxa)
  #calculate summary stats
  L1 <- median(trace_1$Likelihood)
  L2 <- median(trace_2$Likelihood)
  k1 <- 21
  k2 <- 10
  
  df$chisse_lnL[i] <- L1
  df$csse_lnL[i] <- L2
  
  df$chisse_AIC[i] <- 2 * k1 - 2 * L1 
  df$csse_AIC[i] <- 2 * k2 - 2 * L2

  df$chisse_AICc[i] <- (2 * k1 - 2 * L1) + (2 * k1 * (k1 + 1)) / (n_taxa - k1 - 1)
  df$csse_AICc[i] <- (2 * k2 - 2 * L2) + (2 * k2 * (k2 + 1)) / (n_taxa - k2 - 1)
  
  df$chisse_BIC[i] <- log(n_taxa) * k1 - 2 * L1
  df$csse_BIC[i] <- log((n_taxa)) * k2 - 2 * L2
}
# calculate dAIC
df[, c(5,7)] <- t(apply(df[, c(4,6)], 1, function(x) x - min(x)))
# calculate dAICc
df[, c(9,11)] <- t(apply(df[, c(8,10)], 1, function(x) x - min(x)))
# calculate dBIC
df[, c(13,15)] <- t(apply(df[, c(12,14)], 1, function(x) x - min(x)))

write.csv(df[complete.cases(df),], file = "/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/3.Results/sims/modelComparisonsTable_simcsse.csv",
          row.names = F)
# write.csv(df, file = "/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/3.Results/sims/modelComparisonsTable_simchisse.csv",
#           row.names = F)


# make model comparison table
setwd("/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/3.Results/sims")
chisse <- read.csv("modelComparisonsTable_simchisse.csv")
chisse <- chisse[,c(5,9,13,7,11,15)]
csse <- read.csv("modelComparisonsTable_simcsse.csv")
csse <- csse[,c(7,11,15,5,9,13)]

table <- data.frame(matrix(nrow = 4, ncol = 8))
colnames(table) <- c("Simulation Model", "Inference Model", "AIC_0", "AICc_0", "BIC_0", "AIC_10", "AICc_10", "BIC_10")
table[,1] <- c("ChromoHiSSE","ChromoHiSSE","ChromoSSE","ChromoSSE")  
table[,2] <- c("ChromoHiSSE", "ChromoSSE", "ChromoHiSSE","ChromoSSE")


# for (i in 1:3) { # what proportion of the time does the data generating model the the lowest AIC, AICc, BIC?
#   
#   chisse_chisse <- mean(chisse[,i] == 0)
#   csse_csse <- mean(csse[,i] == 0)
#   chisse_csse <- 1 - chisse_chisse
#   csse_chisse <- 1 - csse_csse
#   
#   table[,i+2] <- c(chisse_chisse, chisse_csse, csse_chisse, csse_csse)
# }

for (i in 1:6) { # what proportion of the time does the non-data generating model have no support via AIC, AICc, BIC (i.e., dAIC >10)?
if (i <= 3) {
  chisse_chisse <- mean(chisse[,i] == 0)
  csse_csse <- mean(csse[,i] == 0)
  chisse_csse <- 1 - chisse_chisse
  csse_chisse <- 1 - csse_csse
  table[,i+2] <- c(chisse_chisse, chisse_csse, csse_chisse, csse_csse)
  
} else if (i > 3) {
  chisse_chisse <- mean(chisse[,i] >= 10)
  csse_csse <- mean(csse[,i] >= 10)
  chisse_csse <- 1 - chisse_chisse
  csse_chisse <- 1 - csse_csse
  table[,i+2] <- c(chisse_chisse, chisse_csse, csse_chisse, csse_csse)
  
}
}

write.csv(table, "simulation_modelSelectionFrequency.csv")

