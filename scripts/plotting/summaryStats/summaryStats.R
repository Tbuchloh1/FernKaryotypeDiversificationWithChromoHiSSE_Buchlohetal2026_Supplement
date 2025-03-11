library(RevGadgets)
library(tidyverse)
library(HDInterval)
library(common)

workingdir <- "/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/3.Results/chromohisse_v3/output"

## set file paths
# # POLYPODIINEAE
logfile_paths_po <- c("polyp/polyp_chromohisse_tb_1.log",
                   "polyp/polyp_chromohisse_tb_2.log",
                   "polyp/polyp_chromohisse_tb_3.log",
                   "polyp/polyp_chromohisse_tb_4.log",
                   "polyp/polyp_chromohisse_tb_5.log")
# # ASPLENIINEAE
logfile_paths_as <- c("asple/asple_chromohisse_tb_1.log",
                   "asple/asple_chromohisse_tb_2.log",
                   "asple/asple_chromohisse_tb_3.log",
                   "asple/asple_chromohisse_tb_4.log",
                   "asple/asple_chromohisse_tb_5.log")
# # PTERIDINEAE
logfile_paths_pt <- c("pteri/pteri_chromohisse_tb_1.log",
                   "pteri/pteri_chromohisse_tb_2.log",
                   "pteri/pteri_chromohisse_tb_3.log",
                   "pteri/pteri_chromohisse_tb_4.log",
                   "pteri/pteri_chromohisse_tb_5.log")

# Set working directory
setwd(workingdir)

# process data
traces_as <- readTrace(paths = logfile_paths_as, burnin = 0.1)
traces_po <- readTrace(paths = logfile_paths_po, burnin = 0.1)
traces_pt <- readTrace(paths = logfile_paths_pt, burnin = 0.1)

traces_as <- combineTraces(traces_as)
# df_as <- traces_as[[1]]
traces_po <- combineTraces(traces_po)
# df_po <- traces_po[[1]]
traces_pt <- combineTraces(traces_pt)
# df_to <- traces_pt[[1]]

## Ask questions

traces <- list(traces_as, traces_po, traces_pt)
clds <- c("Aspleniineae", "Polypodiineae", "Pteridineae")

# q_lab <- c("Net Diversification, i", "Net Diversification, ii",
#            "T_i","T_ii",
#            "P_i","P_ii",
#            "Polyploid rate difference, i - ii",
#            "Polyploid Cladogenetic, i?","Polyploid Cladogenetic, ii?")

subscript_i <- intToUtf8(0x1D62)
subscript_p <- subsc("p")
q_lab <- c(paste0("r",subscript_i), 
           paste0("r",subscript_i, subscript_i),
           paste0("T",subscript_i),
           paste0("T",subscript_i,subscript_i),
           paste0("P",subscript_i),
           paste0("P",subscript_i,subscript_i),
           paste0("D", subsc("t")),
           paste0("D", subsc("p")),
           paste0("C",subscript_i),
           paste0("C",subscript_i,subscript_i))

rho <- intToUtf8(0x2374)
eta <- intToUtf8(0x03B7)
delta <- intToUtf8(0x03B4)
gamma <- intToUtf8(0x0263)
phi <- intToUtf8(0x03C6)
chi <- intToUtf8(0x03C7)
q_var <- c("netdiv_i","netdiv_ii","T_i","T_ii","P_i","P_ii","total_diffs_i_ii","poly_diffs_i_ii","princ_clado_poly_i","princ_clado_poly_ii")
descriptions <- c("Net Diversification","Net Diversification",
                  paste0("(",gamma,supsc("c"),"+",delta,supsc("c"),"+",rho,supsc("c"),"+",eta,supsc("c"),") - (",chi,supsc("c"),"+",phi,")"),
                  paste0("(",gamma,supsc("c"),"+",delta,supsc("c"),"+",rho,supsc("c"),"+",eta,supsc("c"),") - (",chi,supsc("c"),"+",phi,")"),
                  paste0("(",rho,supsc("c"),"+",eta,supsc("c"),") - (",gamma,supsc("c"),"+",delta,supsc("c"),"+",chi,supsc("c"),"+",phi,")"),
                  paste0("(",rho,supsc("c"),"+",eta,supsc("c"),") - (",gamma,supsc("c"),"+",delta,supsc("c"),"+",chi,supsc("c"),"+",phi,")"),
                  paste0("(",gamma,subscript_i,supsc("c"),"+",delta,subscript_i,supsc("c"),"+",rho,subscript_i,supsc("c"),"+",eta,subscript_i,supsc("c"),") - (",gamma,subscript_i,subscript_i,supsc("c"),"+",delta,subscript_i,subscript_i,supsc("c"),"+",rho,subscript_i,subscript_i,supsc("c"),"+",eta,subscript_i,subscript_i,supsc("c"),")"),
                  paste0("(",rho,subscript_i,supsc("c"),"+",eta,subscript_i,supsc("c"),") - (",rho,subscript_i,subscript_i,supsc("c"),"+",eta,subscript_i,subscript_i,supsc("c"),")"),
                  paste0("(",rho,supsc("c"),"+",eta,supsc("c"),") - (",rho,supsc("a"),"+",eta,supsc("a"),")"),
                  paste0("(",rho,supsc("c"),"+",eta,supsc("c"),") - (",rho,supsc("a"),"+",eta,supsc("a"),")"))

# an empty dataframe for catching outputs
out <- data.frame(matrix(nrow = 1, ncol = 8))
colnames(out) <- c("Clade","Parameter", "Description","Mean","Median","MAP","quantile_2.5","quantile_97.5")

for (i in 1:3) {
  
  ## set current clade's data
  trace <- traces[[i]]
  
  ## calculate summary stats
  # net diversification i
    trace$combined$netdiv_i <- (trace$combined$clado_no_change_a + 
                                  trace$combined$clado_fission_a + 
                                  trace$combined$clado_fusion_a + 
                                  trace$combined$clado_polyploid_a + 
                                  trace$combined$clado_demipoly_a + 
                                  trace$combined$chi) - trace$combined$turnover
  # net diversification ii
    trace$combined$netdiv_ii <- (trace$combined$clado_no_change_b + 
                                   trace$combined$clado_fission_b + 
                                   trace$combined$clado_fusion_b + 
                                   trace$combined$clado_polyploid_b + 
                                   trace$combined$clado_demipoly_b + 
                                   trace$combined$chi) - trace$combined$turnover
    # hist(trace$combined$netdiv_ii)
    # summarizeTrace(trace = trace, vars = "netdiv_ii")
  # Ratio of total lambda from chromosome change in i 
    trace$combined$T_i <- (trace$combined$clado_fission_a + 
                             trace$combined$clado_fusion_a + 
                             trace$combined$clado_polyploid_a + 
                             trace$combined$clado_demipoly_a) - 
                            (trace$combined$clado_no_change_a + 
                               trace$combined$chi)
  # Ratio of total lambda from chromosome change in ii 
    trace$combined$T_ii <- (trace$combined$clado_fission_b + 
                              trace$combined$clado_fusion_b + 
                              trace$combined$clado_polyploid_b + 
                              trace$combined$clado_demipoly_b) - 
                            (trace$combined$clado_no_change_b + 
                               trace$combined$chi)
  
  # Ratio of total lambda from WGD change in i 
    trace$combined$P_i <- (trace$combined$clado_polyploid_a + 
                             trace$combined$clado_demipoly_a) - 
                          (trace$combined$clado_fusion_a + 
                             trace$combined$clado_fission_a + 
                             trace$combined$clado_no_change_a + 
                             trace$combined$chi)
  # Ratio of total lambda from WGD change in ii
    trace$combined$P_ii <- (trace$combined$clado_polyploid_b + 
                              trace$combined$clado_demipoly_b) - 
                           (trace$combined$clado_fusion_b + 
                              trace$combined$clado_fission_b + 
                              trace$combined$clado_no_change_b + 
                              trace$combined$chi)
    
  # Difference between state i and ii in rates of cladogenetic polyploidy
    trace$combined$total_diffs_i_ii <- (trace$combined$clado_fusion_a + 
                                          trace$combined$clado_fission_a + 
                                          trace$combined$clado_polyploid_a + 
                                          trace$combined$clado_demipoly_a) - 
                                        (trace$combined$clado_fusion_b + 
                                           trace$combined$clado_fission_b + 
                                           trace$combined$clado_polyploid_b + 
                                           trace$combined$clado_demipoly_b)
    
  # Difference between state i and ii in rates of cladogenetic polyploidy
    trace$combined$poly_diffs_i_ii <- (trace$combined$clado_polyploid_a + 
                                         trace$combined$clado_demipoly_a) - 
                                      (trace$combined$clado_polyploid_b + 
                                         trace$combined$clado_demipoly_b)
  
  # Test for whether polyploidy is principally cladogenetic in state i
    trace$combined$princ_clado_poly_i <- (trace$combined$clado_polyploid_a + 
                                            trace$combined$clado_demipoly_a) - 
                                         (trace$combined$rho_a + 
                                            trace$combined$eta_a)
  # Test for whether polyploidy is principally cladogenetic in state ii
    trace$combined$princ_clado_poly_ii <- (trace$combined$clado_polyploid_b + 
                                             trace$combined$clado_demipoly_b) - 
                                         (trace$combined$rho_b + 
                                            trace$combined$eta_b)

    
  ## write values to a table
  for (k in 1:length(q_var)) {
    temp_out <- data.frame(summarizeTrace(trace, vars = q_var[k]))
    temp_out <- c(clds[i], q_lab[k], descriptions[k], temp_out[1:5,])
    out <- rbind(out, temp_out)
    out <- out[complete.cases(out),]
  }
}
  write.csv(out, "summaryStatsTable_test.csv", row.names = F, fileEncoding = "UTF-8")

##################
