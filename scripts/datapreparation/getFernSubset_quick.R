# set working directory
setwd("/path/to/working/directory")

############################

# Set relevant file paths
FTOL_tree               <- "polypodiopsida/polypodiopsida_csomepruned.tree"
FTOL_cleanedCounts      <- "polypodiopsida/polypodiopsida_ccdb_mincounts.tsv"
new_tree_file           <- "polypodiidae_exAdigitata/polypodiidae_exAdigitata_csomepruned.tree"
new_counts_file         <- "polypodiidae_exAdigitata/polypodiidae_exAdigitata_counts.tsv"
new_expandedcounts_file <- "polypodiidae_exAdigitata/polypodiidae_exAdigitata_counts_expanded.tsv"

############################

# Provide tip labels determining your targeted clade
tip_label_1 <- "Osmunda_regalis"
tip_label_2 <- "Acrostichum_danaeifolium"

############################

## prune the FTOL to the focal clade(s)
library(phytools)
library(ape)
tree <- read.tree(FTOL_tree) # load the pruned FTOL
# str(tree)

# extract focal node (tips hand-picked by visual inspection of full)
MRCA <- findMRCA(tree, tips = c(tip_label_1,tip_label_2), type = c("node"))

# extract all descendant nodes and tips
desc <- getDescendants(tree, node = MRCA) 

# get the tip labels
tips <- tree$tip.label[desc] 

# keep only the tips in the pruned FTOL tree
tips_comp <- tips[complete.cases(tree$tip.label[desc])] 

# prune the FTOL tree to focal clade
tree_subset <- drop.tip(tree, tree$tip.label[-na.omit(match(tips_comp, tree$tip.label))]) 

# write the pruned tree
write.tree(tree_subset, file = new_tree_file) 

## prune chromosome counts to the focal clades
ccdb.ftol <- read.table(FTOL_cleanedCounts, sep = '\t', header = FALSE) # read all FTOL counts
ccdb.subset <- ccdb.ftol[ccdb.ftol$V1 %in% tree_subset$tip.label,] # prune counts to subtree

## make final data frames for exporting
# change the data frame to a *named vector* to make final .tsv for non-hidden state analyses
counts <- as.numeric(ccdb.subset[,2])

names(counts) <- ccdb.subset[,1]

counts <- data.frame(counts)

# Modify a new .tsv for hidden state analyses
# expand data frame for hidden states
# max chromosomes + 10 more (Mayrose et al. 2010 recommendation)
maximum.count <- max(ccdb.subset[,2]) + 10 

second_column <- as.numeric(ccdb.subset[,2])

expand_data <- second_column+maximum.count+1 # Expansion (leave room for '0' in hidden state )

# creating the double entry in each datum
counts.expanded <- paste0("(", second_column, " ", expand_data,")")

names(counts.expanded) <- ccdb.subset[,1]

counts.expanded <- data.frame(counts.expanded)

# write simple data to .tsv
write.table(counts, 
            file = new_counts_file,
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)

# write hidden state data to .tsv
write.table(counts.expanded, 
            file = new_expandedcounts_file,
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)

