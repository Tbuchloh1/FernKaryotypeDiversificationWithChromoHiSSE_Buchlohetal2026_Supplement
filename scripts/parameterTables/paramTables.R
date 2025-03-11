library(RevGadgets)
library(coda)

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
traces1 <- combineTraces(traces)
# df1 <- traces[[1]]

# ASPLENIINEAE
logfile_paths <- c("chromohisse_v3/output/asple/asple_chromohisse_tb_1.log",
                   "chromohisse_v3/output/asple/asple_chromohisse_tb_2.log",
                   "chromohisse_v3/output/asple/asple_chromohisse_tb_3.log",
                   "chromohisse_v3/output/asple/asple_chromohisse_tb_4.log",
                   "chromohisse_v3/output/asple/asple_chromohisse_tb_5.log")
# process data
traces <- readTrace(paths = logfile_paths, burnin = 0.1)
traces2 <- combineTraces(traces)
# df2 <- traces[[1]]

# PTERIDINEAE
logfile_paths <- c("chromohisse_v3/output/pteri/pteri_chromohisse_tb_1.log",
                   "chromohisse_v3/output/pteri/pteri_chromohisse_tb_2.log",
                   "chromohisse_v3/output/pteri/pteri_chromohisse_tb_3.log",
                   "chromohisse_v3/output/pteri/pteri_chromohisse_tb_4.log",
                   "chromohisse_v3/output/pteri/pteri_chromohisse_tb_5.log")
# process data
traces <- readTrace(paths = logfile_paths, burnin = 0.1)
traces3 <- combineTraces(traces)
# df3 <- traces[[1]]


# # CHROMOSSE
# POLYPODIINEAE
logfile_paths <- c("chromosse_v3/output/polyp/polyp_chromosse_tb_1.log",
                   "chromosse_v3/output/polyp/polyp_chromosse_tb_2.log",
                   "chromosse_v3/output/polyp/polyp_chromosse_tb_3.log",
                   "chromosse_v3/output/polyp/polyp_chromosse_tb_4.log",
                   "chromosse_v3/output/polyp/polyp_chromosse_tb_5.log")
# process data
traces <- readTrace(paths = logfile_paths, burnin = 0.1)
traces4 <- combineTraces(traces)
# df4 <- traces[[1]]

# ASPLENIINEAE
logfile_paths <- c("chromosse_v3/output/asple/asple_chromosse_tb_1.log",
                   "chromosse_v3/output/asple/asple_chromosse_tb_2.log",
                   "chromosse_v3/output/asple/asple_chromosse_tb_3.log",
                   "chromosse_v3/output/asple/asple_chromosse_tb_4.log",
                   "chromosse_v3/output/asple/asple_chromosse_tb_5.log")
# process data
traces <- readTrace(paths = logfile_paths, burnin = 0.1)
traces5 <- combineTraces(traces)
# df5 <- traces[[1]]

# # PTERIDINEAE
logfile_paths <- c("chromosse_v3/output/pteri/pteri_chromosse_tb_1.log",
                   "chromosse_v3/output/pteri/pteri_chromosse_tb_2.log",
                   "chromosse_v3/output/pteri/pteri_chromosse_tb_3.log",
                   "chromosse_v3/output/pteri/pteri_chromosse_tb_4.log",
                   "chromosse_v3/output/pteri/pteri_chromosse_tb_5.log")
# process data
traces <- readTrace(paths = logfile_paths, burnin = 0.1)
traces6 <- combineTraces(traces)
# df6 <- traces[[1]]


tracesN <- list(traces1,traces2,traces3,traces4,traces5,traces6)
vars_ch <- c("chi",
          "clado_fission_a","clado_fission_b",
          "clado_fusion_a","clado_fusion_b",
          "clado_polyploid_a","clado_polyploid_b",
          "clado_demipoly_a","clado_demipoly_b",
          "clado_no_change_a","clado_no_change_b",
          "alpha",
          "gamma_a","gamma_b",
          "delta_a","delta_b",
          "rho_a","rho_b",
          "eta_a","eta_b",
          "turnover")
vars_c <- c("clado_fission",
            "clado_fusion",
            "clado_polyploid",
            "clado_demipoly",
            "clado_no_change",
            "gamma",
            "delta",
            "rho",
            "eta",
            "turnover")
clades <- rep(c("Polyp","Asple","Pteri"),2)
for (i in 1:6) {
  traces_temp <- tracesN[[i]]
  clade <- clades[i]
  
  if (i <= 3) {
  j <- summarizeTrace(traces_temp,vars=vars_ch)
  df <- data.frame(matrix(nrow = length(vars_ch), ncol = length(j$chi$combined)+1))
  colnames(df) <- c("parameter", names(j$chi$combined))
  df$parameter <- vars_ch
  
  for (k in 1:length(vars_ch)) {
    temp <- t(as.data.frame(j[[k]]))
    df[k,2:ncol(df)] <- temp[1,]
  }
  filename = paste0("chromohisse_v3/outputTables/",clade,"_chromohisse_parameterestimates.csv")
  write.csv(df, file = filename)
  } else if (i > 3) {
    j <- summarizeTrace(traces_temp,vars=vars_c)
    df <- data.frame(matrix(nrow = length(vars_c), ncol = length(j$rho$combined)+1))
    colnames(df) <- c("parameter", names(j$rho$combined))
    df$parameter <- vars_c
    
    for (k in 1:length(vars_c)) {
      temp <- t(as.data.frame(j[[k]]))
      df[k,2:ncol(df)] <- temp[1,]
    }
    filename = paste0("chromosse_v3/outputTables/",clade,"_chromosse_parameterestimates.csv")
    write.csv(df, file = filename)
  }
} 
 

