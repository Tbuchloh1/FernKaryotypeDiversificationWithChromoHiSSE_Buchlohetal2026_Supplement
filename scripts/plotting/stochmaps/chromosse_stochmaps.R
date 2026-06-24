# library(devtools)
# devtools::install_github("revbayes/RevGadgets@stochastic_map")
library(RevGadgets) #dependency yulab.utils not updated when installing revgadgets with dependencies. 
library(ape)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(treeio)

workingdir <- "/path/to/working/directory"

# POLYPODIINEAE
# tree_path <- "1.Data/polypodiineae/polypodiineae_csomepruned.tree"
# counts <- "1.Data/polypodiineae/polypodiineae_counts.tsv"
# stochMaps <- c("3.Results/chromosse_v3/output/polyp/polyp_chromosse_tb_SCM_1.log",
#                "3.Results/chromosse_v3/output/polyp/polyp_chromosse_tb_SCM_2.log",
#                "3.Results/chromosse_v3/output/polyp/polyp_chromosse_tb_SCM_3.log",
#                "3.Results/chromosse_v3/output/polyp/polyp_chromosse_tb_SCM_4.log",
#                "3.Results/chromosse_v3/output/polyp/polyp_chromosse_tb_SCM_5.log")
# max_k = 174


# # PTERIDINEAE
# tree_path <- "1.Data/pteridaceae/pteridaceae_csomepruned.tree"
# counts <- "1.Data/pteridaceae/pteridaceae_counts.tsv"
# stochMaps <- c("3.Results/chromosse_v3/output/pteri/pteri_chromosse_tb_SCM_1.log",
#                "3.Results/chromosse_v3/output/pteri/pteri_chromosse_tb_SCM_2.log",
#                "3.Results/chromosse_v3/output/pteri/pteri_chromosse_tb_SCM_3.log",
#                "3.Results/chromosse_v3/output/pteri/pteri_chromosse_tb_SCM_4.log",
#                "3.Results/chromosse_v3/output/pteri/pteri_chromosse_tb_SCM_5.log")
# max_k = 238

# ASPLENIINEAE
tree_path <- "1.Data/aspleniineae/aspleniineae_csomepruned.tree"
counts <- "1.Data/aspleniineae/aspleniineae_counts.tsv"
stochMaps <- c("3.Results/chromosse_v3/output/asple/asple_chromosse_tb_SCM_1.log",
               "3.Results/chromosse_v3/output/asple/asple_chromosse_tb_SCM_2.log",
               "3.Results/chromosse_v3/output/asple/asple_chromosse_tb_SCM_3.log",
               "3.Results/chromosse_v3/output/asple/asple_chromosse_tb_SCM_4.log",
               "3.Results/chromosse_v3/output/asple/asple_chromosse_tb_SCM_5.log")
max_k = 174

# set working directory
setwd(workingdir)

##### get data ##### 
# read in tree 
tree <- RevGadgets::readTrees(tree_path)
cts <- max_k
# read in and process maps
maps <- processStochMaps(tree = tree,
                         paths = stochMaps,
                         states = as.character(1:cts)) #double if using hidden state model

# # rename states #not necessary for non-hidden model
# cns <- c(c("node", "bl", "x0", "x1", "vert"),
#          paste0(1:cts, "A"),
#          paste0(1:cts, "B"))
# colnames(maps) <- cns

# set up tip data
chromo_data <- read.table(counts, sep = "\t") #read in table of counts
chromo_data_vec <- chromo_data$V2 # make vector of chromosome counts
names(chromo_data_vec) <- chromo_data$V1 # make it a named vector of counts

# ladderize the tip data so it plots right
tree_temp <- read.tree(tree_path)
tree_temp <- as.phylo(tree_temp)
tree_temp <- ladderize(tree_temp, right = F)
ladderized_tip_order <- tree_temp$edge[tree_temp$edge[, 2] <= length(tree_temp$tip.label), 2]
ladderized_tip_labels <- tree_temp$tip.label[ladderized_tip_order]
dat <- chromo_data_vec[ladderized_tip_labels]
# dat <- chromo_data_vec[tree[[1]][[1]]@phylo$tip.label] # reorder the named vector to match the order of the tree object. BROKEN
dat_rel <- dat / max(dat) # relative number

tip_dat <- data.frame(dat = dat + 5,
                      dat_rel = dat_rel,
                      order = 1:nrow(chromo_data))

# split into AB maps and chromosome number maps

# # AB
# maps_AB <- maps[ , 1:5]
# maps_AB$A <- rowSums(maps[ , paste0(1:cts, "A")])
# maps_AB$B <- rowSums(maps[ , paste0(1:cts, "B")])

#chromos
# mat <- matrix(nrow = nrow(maps), ncol = cts)
# colnames(mat) <- as.character(1:cts)
# for (i in 1:cts){
#   columns <- i
#   mat[ , i] <- rowSums(maps[, columns])
# }
# 
# maps_chromos <- cbind(maps[ , 1:5], 
#                       mat)

# # make AB map plot 
# hidden <- plotStochMaps(tree[[1]][[1]],
#                         maps_AB,
#                         tip_labels = F) +
#   theme(legend.position = 'none',
#         legend.title = element_text(size = 25),
#         legend.text = element_text(size = 20),
#         legend.key.width = unit(100, "pt")) + 
#   theme(panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank())
# 
# make chromosome number map plot
colors <- viridis::turbo(cts)
names(colors) <- as.character(1:cts)

chromos <- plotStochMaps(tree,
                         maps,
                         tip_labels = F,
                         colors = colors) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank())

# set up fonts
extrafont::font_import() #do once 
extrafont::loadfonts(device = c("all"))

# make tip data barplot 
tips <- ggplot(tip_dat) +
  geom_bar(aes(x = order, y = dat),
           stat="identity") +
  scale_x_continuous(expand = c(0.01, 0)) +
  ggthemes::theme_few() +
  ylab("# chromosomes (n)") +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 10,
                                    family = "Times New Roman"),
        axis.line.y = element_line(color = "grey30"),
        axis.text.y = element_text(size = 5,
                                   family = "Times New Roman"),
        plot.margin = unit(c(0,0.15,0,0.9), "cm"))

lay <- rbind(c(1),
             c(2),
             c(2))

# pdf("figures/trees/trees.pdf", height = 20, width = 20)
grid.arrange(tips, 
             chromos + ggpubr::rotate() + theme(plot.margin = unit(c(0,10,0,57), "pt")), 
             nrow = 2,
             layout_matrix = lay)
# dev.off()
