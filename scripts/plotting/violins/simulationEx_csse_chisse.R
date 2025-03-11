library(RevGadgets)
library(ggplot2)
library(tidyverse)

workingdir <- "/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/3.Results/sims/sim_reduced/output/csse_chisse"

## set file paths
logfile_paths <- "csse_chisse_110.log"



# Set working directory
setwd(workingdir)

# process data
traces <- readTrace(paths = logfile_paths, burnin = 2000)

# traces <- combineTraces(traces)
df <- traces[[1]]

# make separate ana and clado dfs 
clado_cols <- c("chi", 
                "clado_fission_a", "clado_fission_b",
                "clado_fusion_a", "clado_fusion_b",
                "clado_polyploid_a", "clado_polyploid_b",
                "clado_demipoly_a", "clado_demipoly_b",
                "clado_no_change_a", "clado_no_change_b")
ana_cols <- c("alpha",
              "gamma_a", "gamma_b",
              "delta_a", "delta_b",
              "rho_a", "rho_b",
              "eta_a", "eta_b")

df_clado <- df[,clado_cols]
df_ana <- df[,ana_cols]
# df_turnover <- df$turnover

# make summary rate df 
# df_summ <- data.frame(# A
#   diff_change_A = df_clado$clado_fission_a + df_clado$clado_fusion_a + df_clado$clado_polyploid_a + df_clado$clado_demipoly_a - # speciation and chromo change
#     df_clado$clado_no_change_a - df_clado$chi,# speciation and no chromo change
#   # B
#   diff_change_B = df_clado$clado_fission_b + df_clado$clado_fusion_b + df_clado$clado_polyploid_b + df_clado$clado_demipoly_b - # speciation and chromo change
#     df_clado$clado_no_change_b - df_clado$chi, # speciation and no chromo change
#   
#   speciation_total_A = df_clado$clado_fission_a + df_clado$clado_fusion_a + df_clado$clado_polyploid_a + df_clado$clado_demipoly_a + df_clado$clado_no_change_a + df_clado$chi, # sum all clado rates
#   
#   speciation_total_B = df_clado$clado_fission_b + df_clado$clado_fusion_b + df_clado$clado_polyploid_a + df_clado$clado_demipoly_a + df_clado$clado_no_change_b + df_clado$chi # sum all clado rates
#   )

# # Thomas playing around with ideas
# demipolyploidy_diff_BA = df_clado$clado_demipoly_b - df_clado$clado_demipoly_a
# hist(df_summ$demipolyploidy_diff_BA)
# hdi(df_summ$demipolyploidy_diff_BA, credMass = 0.89)

colnames(df_clado) <- c("Cladogenetic Change in Hidden",
                        "Cladogenetic Fission A",
                        "Cladogenetic Fission B",
                        "Cladogenetic Fusion A",
                        "Cladogenetic Fusion B",
                        "Cladogenetic Polyploidy A",
                        "Cladogenetic Polyploidy B",
                        "Cladogenetic Demipolyploidy A",
                        "Cladogenetic Demipolyploidy B",
                        "No Cladogenetic Change A",
                        "No Cladogenetic Change B")

colnames(df_ana) <- c("Anagenetic Change in Hidden",
                      "Anagenetic Fission A",
                      "Anagenetic Fission B",
                      "Anagenetic Fusion A",
                      "Anagenetic Fusion B",
                      "Anagenetic Polyploidy A",
                      "Anagenetic Polyploidy B",
                      "Anagenetic Demipolyploidy A",
                      "Anagenetic Demipolyploidy B")


# get summaries for table 
sumC <- summarizeTrace(traces, vars = clado_cols)
sumA <- summarizeTrace(traces, vars = ana_cols)
# sumT <- summarizeTrace(traces, vars = "turnover")

df_clado %>%
tidyr::gather(key = "grp", 
              value = "val",
              factor_key = TRUE) -> df_clado
  
df_ana %>%
  tidyr::gather(key = "grp", 
                value = "val",
                factor_key = TRUE) -> df_ana

# df_summ %>%
#   dplyr::select(diff_change_A, diff_change_B) %>%
#   tidyr::gather(key = "grp", 
#                 value = "val",
#                 factor_key = TRUE) -> df_diff
# 
# df_summ %>%
#   dplyr::select(speciation_total_A, speciation_total_B) %>%
#   tidyr::gather(key = "grp", 
#                 value = "val",
#                 factor_key = TRUE) -> df_speciation

df_clado$grp2 <- car::recode(df_clado$grp, " 'Cladogenetic Change in Hidden' = 'Chi_c'; 
                                             c('Cladogenetic Fission A','Cladogenetic Fusion A','Cladogenetic Polyploidy A',
                                             'Cladogenetic Demipolyploidy A','No Cladogenetic Change A') = 'A'; 
                                             c('Cladogenetic Fission B', 'Cladogenetic Fusion B','Cladogenetic Polyploidy B',
                                             'Cladogenetic Demipolyploidy B','No Cladogenetic Change B') = 'B' ")

df_ana$grp2 <- car::recode(df_ana$grp, " 'Anagenetic Change in Hidden' = 'Chi_a'; 
                                          c('Anagenetic Fission A','Anagenetic Fusion A',
                                          'Anagenetic Polyploidy A','Anagenetic Demipolyploidy A') = 'A'; 
                                          c('Anagenetic Fission B','Anagenetic Fusion B',
                                          'Anagenetic Polyploidy B','Anagenetic Demipolyploidy B') = 'B' ")

##### Plot #####
# all_colors <- c("#9661A8", "#1AE4B6FF", "#FABA39FF", "#CB2A04FF")
all_colors <- c("#0cb2af", "#a1c65d", "#fac723", "#f29222", "#e95e50", "#936fac")
# set up colors
colors <- c(all_colors[1],
            all_colors[2], all_colors[2], 
            all_colors[3], all_colors[3],  
            all_colors[4], all_colors[4],
            all_colors[5], all_colors[5],  
            all_colors[6], all_colors[6])

names(colors) <- levels(df_clado$grp)

################# Specifying sim values ################
manual_vals <- c(NA,0.0017,NA,0.012,NA,0.0582,NA,0.017,NA,0.537,NA)
manual_levels <- levels(df_clado$grp)
manual_df <- data.frame(manual_levels,manual_vals)
manual_df$manual_levels <- as.factor(manual_df$manual_levels)
################# Specifying sim values ################

g <- ggplot(df_clado) +
  geom_hline(yintercept = 0.0, color = "grey") + 
  geom_violin(data = df_clado,
              aes(x = grp,
                  y = val,
                  group = grp,
                  fill = grp,
                  linetype = grp2),
              color = "black",
              lwd = 0.5,
              scale = "width",
              show.legend = F) +
  stat_summary(fun=mean, aes(x = grp, y = val), geom="point", size=2, color="black") +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = c("solid","dashed","blank")) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(name = element_blank(),
                   labels = c(
                              expression(chi^"c"),
                              expression(gamma["i"]^"c"),
                              expression(gamma["ii"]^"c"),
                              expression(delta["i"]^"c"),
                              expression(delta["ii"]^"c"),
                              expression(rho["i"]^"c"),
                              expression(rho["ii"]^"c"),
                              expression(eta["i"]^"c"),
                              expression(eta["ii"]^"c"),
                              expression(phi["i"]^"c"),
                              expression(phi["ii"]^"c"))) +
  ylab("Parameter value") +
  ggthemes::theme_few() +
  theme(axis.text.x = element_text(face="bold", 
                                   size=15, 
                                   hjust = .5),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(face = "bold", 
                                    size = 20),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_point(data = manual_df, aes(x = manual_levels, y = manual_vals), color = 'red', size = 1)

colors <- c(all_colors[1],
            all_colors[2], all_colors[2], 
            all_colors[3], all_colors[3],  
            all_colors[4], all_colors[4],
            all_colors[5], all_colors[5])

################# Specifying sim values ################
manual_vals <- c(NA,0.001,NA,0.0183,NA,0.00313,NA,0.00228,NA)
manual_levels <- levels(df_ana$grp)
manual_df <- data.frame(manual_levels,manual_vals)
manual_df$manual_levels <- as.factor(manual_df$manual_levels)
################# Specifying sim values ################

h <- ggplot(df_ana) +
  geom_hline(yintercept = 0.0, color = "grey") + 
  geom_violin(data = df_ana,
              aes(x = grp, 
                  y = val, 
                  group = grp, 
                  fill = grp,
                  linetype = grp2),
              color = "black",
              lwd = 0.5,
              scale = "width",
              show.legend = F) +
  stat_summary(fun=mean, aes(x = grp, y = val), geom="point", size=2, color="black") +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = c("solid","dashed","blank")) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(name = element_blank(),
                   labels = c(
                              expression(chi^"a"),
                              expression(gamma["i"]^"a"),
                              expression(gamma["ii"]^"a"),
                              expression(delta["i"]^"a"),
                              expression(delta["ii"]^"a"),
                              expression(rho["i"]^"a"),
                              expression(rho["ii"]^"a"),
                              expression(eta["i"]^"a"),
                              expression(eta["ii"]^"a"))) +
  ylab("Parameter value") +
  ggthemes::theme_few() +
  theme(axis.text.x = element_text(face="bold", 
                                   size=15, 
                                   hjust = .5),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(face = "bold", 
                                    size = 20),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))  +
  geom_point(data = manual_df, aes(x = manual_levels, y = manual_vals), color = 'red', size = 1)


lay <- rbind(c(1,1,1,1,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,2,2,2,2,2,2,2))

# pdf("figures/rates/rates.pdf", width = 10, height = 7)
gridExtra::grid.arrange(g, h,
                        layout_matrix = lay)
# dev.off()


