library(RevGadgets)
library(tidyverse)

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
# # df_as <- traces_as[[1]]
traces_po <- combineTraces(traces_po)
# # df_po <- traces_po[[1]]
traces_pt <- combineTraces(traces_pt)
# # df_to <- traces_pt[[1]]
traces <- list(traces_as, traces_po, traces_pt)
clades <- c("Polypodiineae","Aspleniineae","Pteridineae")
## Ask questions

df <- data.frame(matrix(nrow = 0, ncol = 4))


for (i in 1:3) {
  
  ## set current clade's data
  trace <- traces[[i]]
  trace <- trace[[1]]
  n <- nrow(trace)
  temp <- data.frame(matrix(nrow = 0, ncol = 4))
  ## calculate summary stats
  # net diversification i
    temp[1:n,1] <- rep(clades[i], n)
    temp[1:n,2] <- rep("Net Diversification (r)", n)
    temp[1:n,3] <- rep("i", n)
    temp[1:n,4] <- (trace$clado_no_change_a + 
                 trace$clado_fission_a + 
                 trace$clado_fusion_a + 
                 trace$clado_polyploid_a + 
                 trace$clado_demipoly_a + 
                 trace$chi) - trace$turnover
   df <- rbind(df, temp)
  # net diversification ii
    temp[1:n,1] <- rep(clades[i], n)
    temp[1:n,2] <- rep("Net Diversification (r)", n)
    temp[1:n,3] <- rep("ii", n)
    temp[1:n,4] <- (trace$clado_no_change_b + 
                                   trace$clado_fission_b + 
                                   trace$clado_fusion_b + 
                                   trace$clado_polyploid_b + 
                                   trace$clado_demipoly_b + 
                                   trace$chi) - trace$turnover
    df <- rbind(df, temp)
    # hist(trace$combined$netdiv_ii)
    # summarizeTrace(trace = trace, vars = "netdiv_ii")
    
  # Ratio of total lambda from chromosome change in i 
    temp[1:n,1] <- rep(clades[i], n)
    # temp[1:n,2] <- rep("[(ɣᶜ+δᶜ+ρᶜ+ηᶜ) - (χᶜ+φ)] (T)", n)
    temp[1:n,2] <- rep("Relative Karyotype Speciation (RKS)", n)
    temp[1:n,3] <- rep("i", n)
    temp[1:n,4] <- (trace$clado_fission_a + 
                             trace$clado_fusion_a + 
                             trace$clado_polyploid_a + 
                             trace$clado_demipoly_a) - 
                            (trace$clado_no_change_a + 
                               trace$chi)
    df <- rbind(df, temp)
  # Ratio of total lambda from chromosome change in ii 
    temp[1:n,1] <- rep(clades[i], n)
    # temp[1:n,2] <- rep("[(ɣᶜ+δᶜ+ρᶜ+ηᶜ) - (χᶜ+φ)] (T)", n)
    temp[1:n,2] <- rep("Relative Karyotype Speciation (RKS)", n)
    temp[1:n,3] <- rep("ii", n)
    temp[1:n,4] <- (trace$clado_fission_b + 
                              trace$clado_fusion_b + 
                              trace$clado_polyploid_b + 
                              trace$clado_demipoly_b) - 
                            (trace$clado_no_change_b + 
                               trace$chi)
    df <- rbind(df, temp)
  
  # Test for whether polyploidy is principally cladogenetic in state i
    temp[1:n,1] <- rep(clades[i], n)
    # temp[1:n,2] <- rep("[(ρᶜ+ηᶜ) - (ρᵃ+ηᵃ)] (C)", n)
    temp[1:n,2] <- rep("Polyploidy Cladogenetic (PPC)", n)
    temp[1:n,3] <- rep("i", n)
    temp[1:n,4] <- (trace$clado_polyploid_a + trace$clado_demipoly_a) - 
                    (trace$rho_a + trace$eta_a)
    df <- rbind(df, temp)
    
  # Test for whether polyploidy is principally cladogenetic in state ii
    temp[1:n,1] <- rep(clades[i], n)
    # temp[1:n,2] <- rep("[(ρᶜ+ηᶜ) - (ρᵃ+ηᵃ)] (C)", n)
    temp[1:n,2] <- rep("Polyploidy Cladogenetic (PPC)", n)
    temp[1:n,3] <- rep("ii", n)
    temp[1:n,4] <- (trace$clado_polyploid_b + trace$clado_demipoly_b) - 
                    (trace$rho_b + trace$eta_b)
    df <- rbind(df, temp)
        
}

colnames(df) <- c("clade","variable","hid_state","value")

##################

df$variable <- as.factor(df$variable)
df$variable <- factor(df$variable, levels = c("Net Diversification (r)", "Relative Karyotype Speciation (RKS)","Polyploidy Cladogenetic (PPC)"))
df$hid_state <- as.factor(df$hid_state)
# Plot Net div., T, and TCP - TAP

df_filtered <- df %>%
  group_by(clade, variable, hid_state) %>%
  filter(value >= quantile(value, 0.025) & value <= quantile(value, 0.975))

 
ggplot(df_filtered, aes(x = hid_state, y = value, fill = hid_state)) + 
  geom_boxplot() +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  facet_grid(rows = vars(variable), cols = vars(clade), switch = "y", scales = "free_y") +  # Arrange in 3x3 grid
  scale_fill_manual(values = c("i" = "skyblue", "ii" = "salmon")) + 
  theme_classic() +
  theme(legend.position = "none") +  # Remove legend
  labs(x = "Hidden State", y = "Parameter Value") 



  
  
  
  
  
  
  
  
  
  
  
  
