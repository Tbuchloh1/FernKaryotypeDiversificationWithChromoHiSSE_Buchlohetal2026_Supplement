setwd("/path/to/working/directory")

## prune the FTOL to the focal clades: Pteridaceae (222tips), Dryopteridaceae (281tips), Polypodiineae w/o Dryopteridaceae (191tips), and Aspleniineae (267tips), Polypodiidae (1170tips)
library(phytools)
tree <- read.tree("polypodiopsida/polypodiopsida_csomepruned.tree") # load the pruned FTOL
# str(tree)

# extract focal node (tips hand-picked by visual inspection of full)
MRCA_asple <- findMRCA(tree, tips = c("Leptogramma_pozoi","Gymnocarpium_robertianum"), type = c("node"))
MRCA_dryop <- findMRCA(tree, tips = c("Arachniodes_denticulata","Lastreopsis_hispida"), type = c("node")) 
MRCA_lepto <- findMRCA(tree, tips = c("Osmundastrum_cinnamomeum","Acrostichum_danaeifolium"), type = c("node")) 
MRCA_polyp <- findMRCA(tree, tips = c("Adenophorus_hymenophylloides","Lomariopsis_palustris"), type = c("node")) 
MRCA_pteri <- findMRCA(tree, tips = c("Acrostichum_danaeifolium","Pellaea_ovata"), type = c("node")) 
MRCA_polypdry <- findMRCA(tree, tips = c("Didymochlaena_truncatula","Adenophorus_hymenophylloides"), type = c("node")) 

# extract all descendant nodes and tips
desc_asple <- getDescendants(tree, node = MRCA_asple) 
desc_dryop <- getDescendants(tree, node = MRCA_dryop) 
desc_lepto <- getDescendants(tree, node = MRCA_lepto) 
desc_polyp <- getDescendants(tree, node = MRCA_polyp) 
desc_pteri <- getDescendants(tree, node = MRCA_pteri) 
desc_polypdry <- getDescendants(tree, node = MRCA_polypdry) 

# get the tip labels
tips_asple <- tree$tip.label[desc_asple] 
tips_dryop <- tree$tip.label[desc_dryop] 
tips_lepto <- tree$tip.label[desc_lepto] 
tips_polyp <- tree$tip.label[desc_polyp] 
tips_pteri <- tree$tip.label[desc_pteri] 
tips_polypdry <- tree$tip.label[desc_polypdry] 

# keep only the tips in the pruned FTOL tree
tips_comp_asple <- tips_asple[complete.cases(tree$tip.label[desc_asple])] 
tips_comp_dryop <- tips_dryop[complete.cases(tree$tip.label[desc_dryop])] 
tips_comp_lepto <- tips_lepto[complete.cases(tree$tip.label[desc_lepto])] 
tips_comp_polyp <- tips_polyp[complete.cases(tree$tip.label[desc_polyp])] 
tips_comp_pteri <- tips_pteri[complete.cases(tree$tip.label[desc_pteri])] 
tips_comp_polypdry <- tips_polypdry[complete.cases(tree$tip.label[desc_polypdry])] 

# prune the FTOL tree to focal clade
tree_asple <- drop.tip(tree, tree$tip.label[-na.omit(match(tips_comp_asple, tree$tip.label))]) 
tree_dryop <- drop.tip(tree, tree$tip.label[-na.omit(match(tips_comp_dryop, tree$tip.label))]) 
tree_lepto <- drop.tip(tree, tree$tip.label[-na.omit(match(tips_comp_lepto, tree$tip.label))]) 
tree_polyp <- drop.tip(tree, tree$tip.label[-na.omit(match(tips_comp_polyp, tree$tip.label))]) 
tree_pteri <- drop.tip(tree, tree$tip.label[-na.omit(match(tips_comp_pteri, tree$tip.label))]) 
tree_polypdry <- drop.tip(tree, tree$tip.label[-na.omit(match(tips_comp_polypdry, tree$tip.label))]) 

# write the pruned tree
write.tree(tree_asple, file = "aspleniineae/aspleniineae_csomepruned.tree") 
write.tree(tree_dryop, file = "dryopteridaceae/dryopteridaceae_csomepruned.tree") 
write.tree(tree_lepto, file = "polypodiidae/polypodiidae_csomepruned.tree") 
write.tree(tree_polyp, file = "polypodiineae/polypodiineae_csomepruned.tree") 
write.tree(tree_pteri, file = "pteridaceae/pteridaceae_csomepruned.tree") 
write.tree(tree_polypdry, file = "polypodiineae_comp/polypodiineae_comp_csomepruned.tree") 

## prune chromosome counts to the focal clades
ccdb.ftol <- read.table("polypodiopsida/polypodiopsida_ccdb_mincounts.tsv", sep = '\t', header = FALSE) # read all FTOL counts
ccdb.asple <- ccdb.ftol[ccdb.ftol$V1 %in% tree_asple$tip.label,] # prune counts to subtree
ccdb.dryop <- ccdb.ftol[ccdb.ftol$V1 %in% tree_dryop$tip.label,] # prune counts to subtree
ccdb.lepto <- ccdb.ftol[ccdb.ftol$V1 %in% tree_lepto$tip.label,] # prune counts to subtree
ccdb.polyp <- ccdb.ftol[ccdb.ftol$V1 %in% tree_polyp$tip.label,] # prune counts to subtree
ccdb.pteri <- ccdb.ftol[ccdb.ftol$V1 %in% tree_pteri$tip.label,] # prune counts to subtree
ccdb.polypdry <- ccdb.ftol[ccdb.ftol$V1 %in% tree_polypdry$tip.label,] # prune counts to subtree

## make final data frames for exporting
# change the data frame to a *named vector* to make final .tsv for non-hidden state analyses
asple.counts <- as.numeric(ccdb.asple[,2])
dryop.counts <- as.numeric(ccdb.dryop[,2])
lepto.counts <- as.numeric(ccdb.lepto[,2])
polyp.counts <- as.numeric(ccdb.polyp[,2])
pteri.counts <- as.numeric(ccdb.pteri[,2])
polypdry.counts <- as.numeric(ccdb.polypdry[,2])

names(asple.counts) <- ccdb.asple[,1]
names(dryop.counts) <- ccdb.dryop[,1]
names(lepto.counts) <- ccdb.lepto[,1]
names(polyp.counts) <- ccdb.polyp[,1]
names(pteri.counts) <- ccdb.pteri[,1]
names(polypdry.counts) <- ccdb.polypdry[,1]

asple.counts <- data.frame(asple.counts)
dryop.counts <- data.frame(dryop.counts)
lepto.counts <- data.frame(lepto.counts)
polyp.counts <- data.frame(polyp.counts)
pteri.counts <- data.frame(pteri.counts)
polypdry.counts <- data.frame(polypdry.counts)

# Modify a new .tsv for hidden state analyses
# expand data frame for hidden states
# max chromosomes + 10 more (Mayrose et al. 2010 recommendation)
maximum.asple <- max(ccdb.asple[,2]) + 10 
maximum.dryop <- max(ccdb.dryop[,2]) + 10 
maximum.lepto <- max(ccdb.lepto[,2]) + 10 
maximum.polyp <- max(ccdb.polyp[,2]) + 10 
maximum.pteri <- max(ccdb.pteri[,2]) + 10 
maximum.polypdry <- max(ccdb.polypdry[,2]) + 10 

second_column.asple <- as.numeric(ccdb.asple[,2])
second_column.dryop <- as.numeric(ccdb.dryop[,2])
second_column.lepto <- as.numeric(ccdb.lepto[,2])
second_column.polyp <- as.numeric(ccdb.polyp[,2])
second_column.pteri <- as.numeric(ccdb.pteri[,2])
second_column.polypdry <- as.numeric(ccdb.polypdry[,2])

expand_data.asple <- second_column.asple+maximum.asple+1 # Expansion (leave room for '0' in hidden state )
expand_data.dryop <- second_column.dryop+maximum.dryop+1 # Expansion (leave room for '0' in hidden state )
expand_data.lepto <- second_column.lepto+maximum.lepto+1 # Expansion (leave room for '0' in hidden state )
expand_data.polyp <- second_column.polyp+maximum.polyp+1 # Expansion (leave room for '0' in hidden state )
expand_data.pteri <- second_column.pteri+maximum.pteri+1 # Expansion (leave room for '0' in hidden state )
expand_data.polypdry <- second_column.polypdry+maximum.polypdry+1 # Expansion (leave room for '0' in hidden state )

# creating the double entry in each datum
asple.counts.expanded <- paste0("(", second_column.asple, " ", expand_data.asple,")")
dryop.counts.expanded <- paste0("(", second_column.dryop, " ", expand_data.dryop,")")
lepto.counts.expanded <- paste0("(", second_column.lepto, " ", expand_data.lepto,")")
polyp.counts.expanded <- paste0("(", second_column.polyp, " ", expand_data.polyp,")")
pteri.counts.expanded <- paste0("(", second_column.pteri, " ", expand_data.pteri,")")
polypdry.counts.expanded <- paste0("(", second_column.polypdry, " ", expand_data.polypdry,")")

names(asple.counts.expanded) <- ccdb.asple[,1]
names(dryop.counts.expanded) <- ccdb.dryop[,1]
names(lepto.counts.expanded) <- ccdb.lepto[,1]
names(polyp.counts.expanded) <- ccdb.polyp[,1]
names(pteri.counts.expanded) <- ccdb.pteri[,1]
names(polypdry.counts.expanded) <- ccdb.polypdry[,1]

asple.counts.expanded <- data.frame(asple.counts.expanded)
dryop.counts.expanded <- data.frame(dryop.counts.expanded)
lepto.counts.expanded <- data.frame(lepto.counts.expanded)
polyp.counts.expanded <- data.frame(polyp.counts.expanded)
pteri.counts.expanded <- data.frame(pteri.counts.expanded)
polypdry.counts.expanded <- data.frame(polypdry.counts.expanded)

# write simple data to .tsv
write.table(asple.counts, 
            file = "aspleniineae/aspleniineae_counts.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
write.table(dryop.counts, 
            file = "dryopteridaceae/dryopteridaceae_counts.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
write.table(lepto.counts, 
            file = "polypodiidae/polypodiidae_counts.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
write.table(polyp.counts, 
            file = "polypodiineae/polypodiineae_counts.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
write.table(pteri.counts, 
            file = "pteridaceae/pteridaceae_counts.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
write.table(polypdry.counts, 
            file = "polypodiineae_comp/polypodiineae_comp_counts.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)

# write hidden state data to .tsv
write.table(asple.counts.expanded, 
            file = "aspleniineae/aspleniineae_counts_expanded.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
write.table(dryop.counts.expanded, 
            file = "dryopteridaceae/dryopteridaceae_counts_expanded.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
write.table(lepto.counts.expanded, 
            file = "polypodiidae/polypodiidae_counts_expanded.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
write.table(polyp.counts.expanded, 
            file = "polypodiineae/polypodiineae_counts_expanded.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
write.table(pteri.counts.expanded, 
            file = "pteridaceae/pteridaceae_counts_expanded.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
write.table(polypdry.counts.expanded, 
            file = "polypodiineae_comp/polypodiineae_comp_counts_expanded.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
