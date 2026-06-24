library(RevGadgets)
library(ggplot2)
library(tidyverse)

workingdir <- "/path/to/working/directory"

# ASPLENIINEAE
# logfile_paths <- c("asple/asple_chromosse_tb_1.log",
#                    "asple/asple_chromosse_tb_2.log",
#                    "asple/asple_chromosse_tb_3.log",
#                    "asple/asple_chromosse_tb_4.log",
#                    "asple/asple_chromosse_tb_5.log")
# # plot title
# title1 <- "ChromoSSE: Aspleniineae"
# title2 <- "ChromoSSE: Aspleniineae (no phi)"

# POLYPODIINEAE
logfile_paths <- c("polyp/polyp_chromosse_tb_1.log",
                   "polyp/polyp_chromosse_tb_2.log",
                   "polyp/polyp_chromosse_tb_3.log",
                   "polyp/polyp_chromosse_tb_4.log",
                   "polyp/polyp_chromosse_tb_5.log")
# # plot title
title1 <- "ChromoSSE: Polypodiineae"
title2 <- "ChromoSSE: Polypodiineae (no phi)"

# # PTERIDINEAE
# logfile_paths <- c("pteri/pteri_chromosse_tb_1.log",
#                    "pteri/pteri_chromosse_tb_2.log",
#                    "pteri/pteri_chromosse_tb_3.log",
#                    "pteri/pteri_chromosse_tb_4.log",
#                    "pteri/pteri_chromosse_tb_5.log")
# # plot title
# title1 <- "ChromoSSE: Pteridineae"
# title2 <- "ChromoSSE: Pteridineae (no phi)"

setwd(workingdir)

###### ChromoSSE: all parameters ######

# process data
traces <- readTrace(paths = logfile_paths, burnin=0.1)

traces <- combineTraces(traces)
df <- traces[[1]]
cols <- c("gamma","delta", "eta", "rho", "clado_fission", "clado_fusion", "clado_demipoly", "clado_polyploid", "clado_no_change")
#plotTrace(carex, vars = cols)

df <- df[,cols]

df %>%
tidyr::gather(key = "grp", 
              value = "val",
              factor_key = TRUE) -> df


# plot
                 # fission    # fusion   # no change #polyploidy
all_colors <- c("#0cb2af", "#a1c65d", "#fac723", "#f29222", "#936fac")

# set up colors
colors <- c(all_colors[1],
            all_colors[2],
            all_colors[3],
            all_colors[4],
            all_colors[1],
            all_colors[2],
            all_colors[3],
            all_colors[4],
            all_colors[5])

names(colors) <- levels(df$grp)

g <- ggplot(df) +
  geom_hline(yintercept = 0.0, color = "black") + 
  geom_violin(data = df,
              aes(x = grp, 
                  y = val, 
                  group = grp, 
                  fill = grp),
              color = "black",
              lwd = 0.5,
              scale = "width",
              show.legend = F) +
  stat_summary(fun=mean, aes(x = grp, y = val), geom="point", size=2, color="black") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(name = element_blank(),
                   labels = c(expression(gamma^"a"),
                              expression(delta^"a"),
                              expression(eta^"a"),
                              expression(rho^"a"),
                              expression(gamma^"c"),
                              expression(delta^"c"),
                              expression(eta^"c"),
                              expression(rho^"c"),
                              expression(phi^"c"))) +
  ggtitle(title1) +
  ylab("Rate") +
  ggthemes::theme_few() +
  theme(axis.text.x = element_text(face="bold", 
                                   size=20, 
                                   hjust = .5),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(face = "bold", 
                                  size = 25),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# pdf("figures/supplemental_posterior_plots/supp_posteriors_sidrho.pdf", height = 10, width = 12)
g
# dev.off()


###### ChromoSSE: without plotting phi ######

# process data
traces <- readTrace(paths = logfile_paths)

traces <- combineTraces(traces)
df <- traces[[1]]
cols <- c("gamma","delta", "eta", "rho", "clado_fission", "clado_fusion", "clado_demipoly", "clado_polyploid")
#plotTrace(carex, vars = cols)

df <- df[,cols]

df %>%
  tidyr::gather(key = "grp", 
                value = "val",
                factor_key = TRUE) -> df


# plot
# fission    # fusion   # no change #polyploidy
all_colors <- c("#0cb2af", "#a1c65d", "#fac723", "#f29222")

# set up colors
colors <- c(all_colors[1],
            all_colors[2],
            all_colors[3],
            all_colors[4],
            all_colors[1],
            all_colors[2],
            all_colors[3],
            all_colors[4])

names(colors) <- levels(df$grp)

g <- ggplot(df) +
  geom_hline(yintercept = 0.0, color = "black") + 
  geom_violin(data = df,
              aes(x = grp, 
                  y = val, 
                  group = grp, 
                  fill = grp),
              color = "black",
              lwd = 0.5,
              scale = "width",
              show.legend = F) +
  stat_summary(fun=mean, aes(x = grp, y = val), geom="point", size=2, color="black") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(name = element_blank(),
                   labels = c(expression(gamma^"a"),
                              expression(delta^"a"),
                              expression(eta^"a"),
                              expression(rho^"a"),
                              expression(gamma^"c"),
                              expression(delta^"c"),
                              expression(eta^"c"),
                              expression(rho^"c"))) +
  ggtitle(title2) +
  ylab("Rate") +
  ggthemes::theme_few() +
  theme(axis.text.x = element_text(face="bold", 
                                   size=20, 
                                   hjust = .5),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(face = "bold", 
                                  size = 25),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# pdf("figures/supplemental_posterior_plots/supp_posteriors_sidrho.pdf", height = 10, width = 12)
g
# dev.off()
