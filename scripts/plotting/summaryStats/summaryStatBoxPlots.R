library(RevGadgets)
library(tidyverse)
library(stringr)

workingdir <- "/path/to/working/directory"

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
traces_po <- readTrace(paths = logfile_paths_po, burnin = 0.1)
traces_as <- readTrace(paths = logfile_paths_as, burnin = 0.1)
traces_pt <- readTrace(paths = logfile_paths_pt, burnin = 0.1)

traces_po <- combineTraces(traces_po)
# # df_po <- traces_po[[1]]
traces_as <- combineTraces(traces_as)
# # df_as <- traces_as[[1]]
traces_pt <- combineTraces(traces_pt)
# # df_to <- traces_pt[[1]]
traces <- list(traces_po, traces_as, traces_pt)
clades <- c("Polypodiineae","Aspleniineae","Pteridineae")


## Make summaries and plots

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
    temp[1:n,2] <- rep("Karyotype Fraction of Speciation (KFS)", n)
    temp[1:n,3] <- rep("i", n)
    temp[1:n,4] <- (trace$clado_fission_a + 
                             trace$clado_fusion_a + 
                             trace$clado_polyploid_a + 
                             trace$clado_demipoly_a) / 
                    (trace$clado_fission_a + 
                             trace$clado_fusion_a + 
                             trace$clado_polyploid_a + 
                             trace$clado_demipoly_a +
                             trace$clado_no_change_a + 
                             trace$chi)
    df <- rbind(df, temp)
  # Ratio of total lambda from chromosome change in ii 
    temp[1:n,1] <- rep(clades[i], n)
    # temp[1:n,2] <- rep("[(ɣᶜ+δᶜ+ρᶜ+ηᶜ) - (χᶜ+φ)] (T)", n)
    temp[1:n,2] <- rep("Karyotype Fraction of Speciation (KFS)", n)
    temp[1:n,3] <- rep("ii", n)
    temp[1:n,4] <- (trace$clado_fission_b + 
                      trace$clado_fusion_b + 
                      trace$clado_polyploid_b + 
                      trace$clado_demipoly_b) / 
                    (trace$clado_fission_b + 
                       trace$clado_fusion_b + 
                       trace$clado_polyploid_b + 
                       trace$clado_demipoly_b +
                       trace$clado_no_change_b + 
                       trace$chi)
    df <- rbind(df, temp)
  
  # Test for whether polyploidy is principally cladogenetic in state i
    temp[1:n,1] <- rep(clades[i], n)
    # temp[1:n,2] <- rep("[(ρᶜ+ηᶜ) - (ρᵃ+ηᵃ)] (C)", n)
    temp[1:n,2] <- rep("Cladogenetic Fraction of Polyploidy (CFP)", n)
    temp[1:n,3] <- rep("i", n)
    temp[1:n,4] <- (trace$clado_polyploid_a + trace$clado_demipoly_a) / 
                    (trace$clado_polyploid_a + trace$clado_demipoly_a + trace$rho_a + trace$eta_a)
    df <- rbind(df, temp)
  # Test for whether polyploidy is principally cladogenetic in state ii
    temp[1:n,1] <- rep(clades[i], n)
    # temp[1:n,2] <- rep("[(ρᶜ+ηᶜ) - (ρᵃ+ηᵃ)] (C)", n)
    temp[1:n,2] <- rep("Cladogenetic Fraction of Polyploidy (CFP)", n)
    temp[1:n,3] <- rep("ii", n)
    temp[1:n,4] <- (trace$clado_polyploid_b + trace$clado_demipoly_b) / 
                    (trace$clado_polyploid_b + trace$clado_demipoly_b + trace$rho_b + trace$eta_b)
    df <- rbind(df, temp)
        
}

colnames(df) <- c("clade","variable","hid_state","value")

##################

df$variable <- as.factor(df$variable)
df$variable <- factor(df$variable, levels = c("Net Diversification (r)", "Karyotype Fraction of Speciation (KFS)","Cladogenetic Fraction of Polyploidy (CFP)"))
df$hid_state <- as.factor(df$hid_state)
df$clade <- as.factor(df$clade)
df$clade <- factor(df$clade, levels = c("Polypodiineae","Aspleniineae","Pteridineae"))
# Plot Net div., T, and TCP - TAP

df_filtered <- df %>%
  group_by(clade, variable, hid_state) %>%
  filter(value >= quantile(value, 0.025) & value <= quantile(value, 0.975))

# Extract the unique levels of the 'variable' column from your data
unique_variables <- unique(df_filtered$variable)

# Create a data frame with yintercept values for only the existing variables
hline_data <- data.frame(
  variable = unique_variables,  # Use only the variables present in df_filtered
  yintercept = c(0, 0.5, 0.5)    # Specify the yintercept values for each variable
)

# Define y-axis limits for each variable (facet row)
y_axis_limits <- data.frame(
  variable = unique_variables,  # Replace with your actual variable names
  ymin = c(-0.1, 0, 0),                  # Minimum y-axis values for each row
  ymax = c(0.1, 1, 1)                     # Maximum y-axis values for each row
)

# Merge the y-axis limits with your main data
df_filtered <- merge(df_filtered, y_axis_limits, by = "variable")

# Update the ggplot code
ggplot(df_filtered, aes(x = hid_state, y = value, fill = hid_state)) + 
  geom_boxplot() +
  geom_hline(data = hline_data, aes(yintercept = yintercept), linetype = "dashed") +  # Add data-driven hline
  facet_grid(rows = vars(variable), cols = vars(clade), switch = "y", scales = "free_y") +  # Free y-axis
  scale_fill_manual(values = c("i" = "skyblue", "ii" = "salmon")) + 
  theme_classic() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.4),
        panel.spacing.y = unit(0.4, "lines")) +  # adjust spacing between rows if desired)
  labs(x = "Hidden State", y = "Parameter Value") +
  coord_cartesian(ylim = NULL) +  # Use coord_cartesian for dynamic y-axis limits
  geom_blank(aes(y = ymin)) +      # Ensure ymin is respected
  geom_blank(aes(y = ymax))        # Ensure ymax is respected
