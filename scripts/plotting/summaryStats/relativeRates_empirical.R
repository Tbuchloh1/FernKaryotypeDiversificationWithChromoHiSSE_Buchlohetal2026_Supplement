library(RevGadgets)
library(tidyverse)
library(HDInterval)
library(common)

workingdir <- "/path/to/working/directoryt"

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


################## CLADOGENESIS ################## 
# Summarize relative rates across clades to compare consistency of cladogenetic modes of evolution.

params <- c( "rel_chi_a",
             "rel_clado_fission_a", "rel_clado_fusion_a", "rel_clado_polyploid_a", "rel_clado_demipoly_a", "rel_clado_no_change_a", 
             "rel_chi_b",
             "rel_clado_fission_b", "rel_clado_fusion_b", "rel_clado_polyploid_b", "rel_clado_demipoly_b", "rel_clado_no_change_b")

# traces <- list(traces_as, traces_po, traces_pt)
clds <- c("Aspleniineae", "Polypodiineae", "Pteridineae")

df <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(df) <- c("clade", "parameter", "mean", "median", "MAP", "quantile_2.5", "quantile_97.5")

# logfile_paths <- c(logfile_paths_as, logfile_paths_po, logfile_paths_pt)

for (i in 1:3) { # for each clade
  trace <- traces[[i]]
  # trace <- readTrace(paths = logfile_paths[i:(i+4)], burnin = 0.1)
  # trace <- combineTraces(trace)
  clade = clds[i]
  
  # calculate the relative rates
  {  # State i, relative cladogenetic rates for each chromosome parameter
    trace$combined$rel_clado_fission_a <- trace$combined$clado_fission_a /
      (trace$combined$clado_fission_a +
         trace$combined$clado_fusion_a +
         trace$combined$clado_polyploid_a +
         trace$combined$clado_demipoly_a +
         trace$combined$chi +
         trace$combined$clado_no_change_a)
    
    trace$combined$rel_clado_fusion_a <- trace$combined$clado_fusion_a /
      (trace$combined$clado_fission_a +
         trace$combined$clado_fusion_a +
         trace$combined$clado_polyploid_a +
         trace$combined$clado_demipoly_a +
         trace$combined$chi +
         trace$combined$clado_no_change_a)
    
    trace$combined$rel_clado_polyploid_a <- trace$combined$clado_polyploid_a /
      (trace$combined$clado_fission_a +
         trace$combined$clado_fusion_a +
         trace$combined$clado_polyploid_a +
         trace$combined$clado_demipoly_a +
         trace$combined$chi +
         trace$combined$clado_no_change_a)
    
    trace$combined$rel_clado_demipoly_a <- trace$combined$clado_demipoly_a /
      (trace$combined$clado_fission_a +
         trace$combined$clado_fusion_a +
         trace$combined$clado_polyploid_a +
         trace$combined$clado_demipoly_a +
         trace$combined$chi +
         trace$combined$clado_no_change_a)
    
    trace$combined$rel_chi_a <- trace$combined$chi /
      (trace$combined$clado_fission_a +
         trace$combined$clado_fusion_a +
         trace$combined$clado_polyploid_a +
         trace$combined$clado_demipoly_a +
         trace$combined$chi +
         trace$combined$clado_no_change_a)
    
    trace$combined$rel_clado_no_change_a <- trace$combined$clado_no_change_a /
      (trace$combined$clado_fission_a +
         trace$combined$clado_fusion_a +
         trace$combined$clado_polyploid_a +
         trace$combined$clado_demipoly_a +
         trace$combined$chi +
         trace$combined$clado_no_change_a)
    
    # State ii, relative cladogenetic rates for each chromosome parameter
    trace$combined$rel_clado_fission_b <- trace$combined$clado_fission_b /
      (trace$combined$clado_fission_b +
         trace$combined$clado_fusion_b +
         trace$combined$clado_polyploid_b +
         trace$combined$clado_demipoly_b +
         trace$combined$chi +
         trace$combined$clado_no_change_b)
    
    trace$combined$rel_clado_fusion_b <- trace$combined$clado_fusion_b /
      (trace$combined$clado_fission_b +
         trace$combined$clado_fusion_b +
         trace$combined$clado_polyploid_b +
         trace$combined$clado_demipoly_b +
         trace$combined$chi +
         trace$combined$clado_no_change_b)
    
    trace$combined$rel_clado_polyploid_b <- trace$combined$clado_polyploid_b /
      (trace$combined$clado_fission_b +
         trace$combined$clado_fusion_b +
         trace$combined$clado_polyploid_b +
         trace$combined$clado_demipoly_b +
         trace$combined$chi +
         trace$combined$clado_no_change_b)
    
    trace$combined$rel_clado_demipoly_b <- trace$combined$clado_demipoly_b /
      (trace$combined$clado_fission_b +
         trace$combined$clado_fusion_b +
         trace$combined$clado_polyploid_b +
         trace$combined$clado_demipoly_b +
         trace$combined$chi +
         trace$combined$clado_no_change_b)
    
    trace$combined$rel_chi_b <- trace$combined$chi /
      (trace$combined$clado_fission_b +
         trace$combined$clado_fusion_b +
         trace$combined$clado_polyploid_b +
         trace$combined$clado_demipoly_b +
         trace$combined$chi +
         trace$combined$clado_no_change_b)
    
    trace$combined$rel_clado_no_change_b <- trace$combined$clado_no_change_b /
      (trace$combined$clado_fission_b +
         trace$combined$clado_fusion_b +
         trace$combined$clado_polyploid_b +
         trace$combined$clado_demipoly_b +
         trace$combined$chi +
         trace$combined$clado_no_change_b) 
  }
  
  for (j in 1:length(params)) { # within simulations for each parameter (so only have to load traces 1x)
    temp <- c()
    for (k in 1:5) { # within sims and parameters for each summary stat
      temp[k] <- summarizeTrace(trace, vars = params[j])[[1]][[1]][[k]]
      # temp <- c(1:5)
    }
    df[(nrow(df)+1),] <- c(clade, params[j], temp)
  }
}


### Plot cladogenetic estimates (relative) for each parameter across clades 
ggplot(df, aes(x = factor(parameter, levels = params), y = as.numeric(median), color = factor(clade))) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin=as.numeric(quantile_2.5), ymax = as.numeric(quantile_97.5)),
                position = position_dodge(width = 0.5), width = 0.2) +
  theme_classic() +
  labs( x = "", y = "Percent of Cladogenetic Change", color = "Clade") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none") +
  geom_vline(xintercept = 5.5, 
             linetype = "dotted", color = "gray")


################## ANAGENESIS ################## 
# Summarize relative rates across clades to compare consistency of anagenetic modes of evolution.

params <- c( "rel_alpha_a",
             "rel_ana_fission_a", "rel_ana_fusion_a", "rel_ana_polyploid_a", "rel_ana_demipoly_a",
             "rel_alpha_b",
             "rel_ana_fission_b", "rel_ana_fusion_b", "rel_ana_polyploid_b", "rel_ana_demipoly_b")

traces <- list(traces_as, traces_po, traces_pt)
clds <- c("Aspleniineae", "Polypodiineae", "Pteridineae")

df <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(df) <- c("clade", "parameter", "mean", "median", "MAP", "quantile_2.5", "quantile_97.5")

# logfile_paths <- c(logfile_paths_as, logfile_paths_po, logfile_paths_pt)

for (i in 1:3) { # for each clade
  trace <- traces[[i]]
  # trace <- readTrace(paths = logfile_paths[i:(i+4)], burnin = 0.1)
  # trace <- combineTraces(trace)
  clade = clds[i]
  
  # calculate the relative rates
  {  # State i, relative cladogenetic rates for each chromosome parameter
    trace$combined$rel_ana_fission_a <- trace$combined$gamma_a /
      (trace$combined$gamma_a +
         trace$combined$delta_a +
         trace$combined$rho_a +
         trace$combined$eta_a +
         trace$combined$alpha)
    
    trace$combined$rel_ana_fusion_a <- trace$combined$delta_a /
      (trace$combined$gamma_a +
         trace$combined$delta_a +
         trace$combined$rho_a +
         trace$combined$eta_a +
         trace$combined$alpha)
    
    trace$combined$rel_ana_polyploid_a <- trace$combined$rho_a /
      (trace$combined$gamma_a +
         trace$combined$delta_a +
         trace$combined$rho_a +
         trace$combined$eta_a +
         trace$combined$alpha)
    
    trace$combined$rel_ana_demipoly_a <- trace$combined$eta_a /
      (trace$combined$gamma_a +
         trace$combined$delta_a +
         trace$combined$rho_a +
         trace$combined$eta_a +
         trace$combined$alpha)
    
    trace$combined$rel_alpha_a <- trace$combined$alpha /
      (trace$combined$gamma_a +
         trace$combined$delta_a +
         trace$combined$rho_a +
         trace$combined$eta_a +
         trace$combined$alpha)
    
    # State ii, relative cladogenetic rates for each chromosome parameter
    trace$combined$rel_ana_fission_b <- trace$combined$gamma_b /
      (trace$combined$gamma_b +
         trace$combined$delta_b +
         trace$combined$rho_b +
         trace$combined$eta_b +
         trace$combined$alpha)
    
    trace$combined$rel_ana_fusion_b <- trace$combined$delta_b /
      (trace$combined$gamma_b +
         trace$combined$delta_b +
         trace$combined$rho_b +
         trace$combined$eta_b +
         trace$combined$alpha)
    
    trace$combined$rel_ana_polyploid_b <- trace$combined$clado_polyploid_b /
      (trace$combined$gamma_b +
         trace$combined$delta_b +
         trace$combined$rho_b +
         trace$combined$eta_b +
         trace$combined$alpha)
    
    trace$combined$rel_ana_demipoly_b <- trace$combined$eta_a /
      (trace$combined$gamma_b +
         trace$combined$delta_b +
         trace$combined$rho_b +
         trace$combined$eta_b +
         trace$combined$alpha)
    
    trace$combined$rel_alpha_b <- trace$combined$alpha /
      (trace$combined$gamma_b +
         trace$combined$delta_b +
         trace$combined$rho_b +
         trace$combined$eta_b +
         trace$combined$alpha)
  }
  
  for (j in 1:length(params)) { # within simulations for each parameter (so only have to load traces 1x)
    temp <- c()
    for (k in 1:5) { # within sims and parameters for each summary stat
      temp[k] <- summarizeTrace(trace, vars = params[j])[[1]][[1]][[k]]
      # temp <- c(1:5)
    }
    df[(nrow(df)+1),] <- c(clade, params[j], temp)
  }
}


### Plot anagenetic estimates (relative) for each parameter across clades 
ggplot(df, aes(x = factor(parameter, levels = params), y = as.numeric(median), color = factor(clade))) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin=as.numeric(quantile_2.5), ymax = as.numeric(quantile_97.5)),
                position = position_dodge(width = 0.5), width = 0.2) +
  theme_classic() +
  labs( x = "Parameter", y = "Percent of Anagenetic Change", color = "Clade") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_vline(xintercept = 5.5, 
             linetype = "dotted", color = "gray")

