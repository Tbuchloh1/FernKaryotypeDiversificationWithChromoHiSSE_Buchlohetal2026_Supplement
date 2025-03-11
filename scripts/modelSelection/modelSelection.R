library(RevGadgets)

workingdir <- "/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/3.Results"
# Set working directory
setwd(workingdir)

## set file paths
# # CHROMOHISSE
# POLYPODIINEAE
logfile_paths <- c("chromohisse_v3/output/polyp/polyp_chromohisse_tb_1.log",
                   "chromohisse_v3/output/polyp/polyp_chromohisse_tb_2.log",
                   "chromohisse_v3/output/polyp/polyp_chromohisse_tb_3.log",
                   "chromohisse_v3/output/polyp/polyp_chromohisse_tb_4.log",
                   "chromohisse_v3/output/polyp/polyp_chromohisse_tb_5.log")
# process data
traces <- readTrace(paths = logfile_paths, burnin = 0.1)
traces <- combineTraces(traces)
df1 <- traces[[1]]

# ASPLENIINEAE
logfile_paths <- c("chromohisse_v3/output/asple/asple_chromohisse_tb_1.log",
                   "chromohisse_v3/output/asple/asple_chromohisse_tb_2.log",
                   "chromohisse_v3/output/asple/asple_chromohisse_tb_3.log",
                   "chromohisse_v3/output/asple/asple_chromohisse_tb_4.log",
                   "chromohisse_v3/output/asple/asple_chromohisse_tb_5.log")
# process data
traces <- readTrace(paths = logfile_paths, burnin = 0.1)
traces <- combineTraces(traces)
df2 <- traces[[1]]

# PTERIDINEAE
logfile_paths <- c("chromohisse_v3/output/pteri/pteri_chromohisse_tb_1.log",
                   "chromohisse_v3/output/pteri/pteri_chromohisse_tb_2.log",
                   "chromohisse_v3/output/pteri/pteri_chromohisse_tb_3.log",
                   "chromohisse_v3/output/pteri/pteri_chromohisse_tb_4.log",
                   "chromohisse_v3/output/pteri/pteri_chromohisse_tb_5.log")
# process data
traces <- readTrace(paths = logfile_paths, burnin = 0.1)
traces <- combineTraces(traces)
df3 <- traces[[1]]


# # CHROMOSSE
# POLYPODIINEAE
logfile_paths <- c("chromosse_v3/output/polyp/polyp_chromosse_tb_1.log",
                   "chromosse_v3/output/polyp/polyp_chromosse_tb_2.log",
                   "chromosse_v3/output/polyp/polyp_chromosse_tb_3.log",
                   "chromosse_v3/output/polyp/polyp_chromosse_tb_4.log",
                   "chromosse_v3/output/polyp/polyp_chromosse_tb_5.log")
# process data
traces <- readTrace(paths = logfile_paths, burnin = 0.1)
traces <- combineTraces(traces)
df4 <- traces[[1]]

# ASPLENIINEAE
logfile_paths <- c("chromosse_v3/output/asple/asple_chromosse_tb_1.log",
                   "chromosse_v3/output/asple/asple_chromosse_tb_2.log",
                   "chromosse_v3/output/asple/asple_chromosse_tb_3.log",
                   "chromosse_v3/output/asple/asple_chromosse_tb_4.log",
                   "chromosse_v3/output/asple/asple_chromosse_tb_5.log")
# process data
traces <- readTrace(paths = logfile_paths, burnin = 0.1)
traces <- combineTraces(traces)
df5 <- traces[[1]]

# # PTERIDINEAE
logfile_paths <- c("chromosse_v3/output/pteri/pteri_chromosse_tb_1.log",
                   "chromosse_v3/output/pteri/pteri_chromosse_tb_2.log",
                   "chromosse_v3/output/pteri/pteri_chromosse_tb_3.log",
                   "chromosse_v3/output/pteri/pteri_chromosse_tb_4.log",
                   "chromosse_v3/output/pteri/pteri_chromosse_tb_5.log")
# process data
traces <- readTrace(paths = logfile_paths, burnin = 0.1)
traces <- combineTraces(traces)
df6 <- traces[[1]]


## Extract MAP estimates for model likelihood
# set up an empty dataframe for all my outputs
df <- data.frame(matrix(data = NA, nrow = 6, ncol = 11))
colnames(df) <- c("clade","n","model","k","lnL","AIC","dAIC","AICc","dAICc","BIC","dBIC")
df$clade <- c(rep(c("Polyp","Asple","Pteri"),2))
df$n <- rep(c(473,267,222),2)
df$model <- c(rep("ChromHiSSE",3),rep("ChromSSE",3))
df$k <- c(rep("21",3),rep("10",3))

dfs <- list(df1,df2,df3,df4,df5,df6)

# ln(L), AIC, AICc, BIC
for (i in 1:6) {
  df_temp <- dfs[[i]]
  L <- getMAP(df_temp$Likelihood)
  k <- as.numeric(df$k[i])
  AIC_val <- 2 * k - 2 * L 
  AICc_val <- 2 * k - 2 * L + (2 * k * (k + 1)) / (df$n[i] - k - 1)
  BIC_val <- log(df$n[i]) * k - 2 * L
  df$lnL[i] <- L
  df$AIC[i] <- AIC_val
  df$AICc[i] <- AICc_val
  df$BIC[i] <- BIC_val
}

# delta_AIC
for (i in 1:3) {
  ch_aic <- df$AIC[i]
  c_aic <- df$AIC[i+3]
  ch_aicc <- df$AICc[i]
  c_aicc <- df$AICc[i+3]
  ch_bic <- df$BIC[i]
  c_bic <- df$BIC[i+3]
  
  # dAIC
  if (ch_aic > c_aic) {
    ch_dAIC <- ch_aic - c_aic 
    c_dAIC <- 0  
  } else if (c_aic > ch_aic) {
    ch_dAIC <- 0
    c_dAIC <- c_aic - ch_aic
  }
  
  # dAICc
  if (ch_aicc > c_aicc) {
    ch_dAICc <- ch_aicc - c_aicc 
    c_dAICc <- 0  
  } else if (c_aicc > ch_aicc) {
    ch_dAICc <- 0
    c_dAICc <- c_aicc - ch_aicc
  }
  
  # dBIC
  if (ch_bic > c_bic) {
    ch_dBIC <- ch_bic - c_bic 
    c_dBIC <- 0  
  } else if (c_bic > ch_bic) {
    ch_dBIC <- 0
    c_dBIC <- c_bic - ch_bic
  }
  
  # print it to the dataframe
  df$dAIC[i] <- ch_dAIC
  df$dAIC[i+3] <- c_dAIC
  
  df$dAICc[i] <- ch_dAICc
  df$dAICc[i+3] <- c_dAICc
  
  df$dBIC[i] <- ch_dBIC
  df$dBIC[i+3] <- c_dBIC
}
df

write.csv(df, file = "modelComparisonsTable.csv")
  
  


