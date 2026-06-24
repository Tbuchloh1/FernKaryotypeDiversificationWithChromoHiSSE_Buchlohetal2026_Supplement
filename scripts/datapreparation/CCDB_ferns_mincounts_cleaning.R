library(stringr)
library(ftolr)
library(dplyr)
library(ape)

setwd("/path/to/working/directory")

#### load raw CCDB data
raw <- read.csv("1.Data/counts_raw/CCDB_2457Pteridophytes_20240325.csv")
# head(raw)

#### cleaning the dataset from the CCDB 
### very basic cleaning: choose a chromosome count value to care about (in this case, the minimum count as it is the most likely to match the diploid cytotype for a taxon)
## inspect odd-ball counts (by finding min counts less than expected if occasionally triploid (e.g., Dryopteris erythrosora: min=41, median=123))
raw[which(raw$minimum < raw$median/3),]

## apply systematic rule to fix odd-ball minimum counts
# Rule: if odd-ball min count... then, default to median count
rows <- c(which(raw$minimum < raw$median/3)) # vector of row numbers to apply the rule to. 
min_cts <- raw # copy raw so we can write over it without changing the original data frame
for (i in 1:length(rows)) {
  min_cts[rows[i], 4] <- min_cts[rows[i],3] # change in value (col4) to median value (col3)
}
# sanity check: all median and minimum values should now be identical in the following rows
# min_cts[rows,]
# get rid of unneeded columns in the dataframe
n_cts <- min_cts[,4:5]


## more complicated cleaning: names occur once in the tree and only as genus_epithet (to mach tip.label structure (e.g., no var. or subsp. or x))
# splitting 'resolved_name' so we can see the problems and repeats we need to deal with
n_cts[,3:6] <- str_split(n_cts$resolved_name, " ", simplify = T) # new columns 'resolved_name' splits
# view all the issues with epithets
# n_cts[which(n_cts$V5 != ""),]
## fix "Ã\u0097" occurrences
# n_cts[which(n_cts$V4 == "Ã\u0097"),] # view odd-balls
n_cts[which(n_cts$V4 == "Ã\u0097"),4] <- n_cts[which(n_cts$V4 == "Ã\u0097"),5] #shift the epithet into the place of the odd string
## fix occurrences of "x"
# n_cts[which(n_cts$V4 == "×"),] # view odd-balls [NOTE: THAT IS NOT A SIMPLE 'x' it is an '×']
n_cts[which(n_cts$V4 == "×"),4] <- n_cts[which(n_cts$V4 == "×"),5] 
## remove var., f., and subsp. (by removing everything in columns 5 and 6)
n_cts_dups <- n_cts[,1:4]

## remove duplicate names, defaulting to the lowest count. 
# create indicator column
n_cts_dups$dups <- c("n")
# label duplicates
n_cts_dups[duplicated(n_cts_dups[,3:4]) | duplicated(n_cts_dups[,3:4], fromLast = T), 5] <- "y"
n_cts_dups$fullname <- paste(n_cts_dups$V3,n_cts_dups$V4, sep = "_")
# subset duplicate and non-duplicate rows
dups <- n_cts_dups[which(n_cts_dups$dups == "y"),]
not_dups <- n_cts_dups[which(n_cts_dups$dups == "n"),]
# loop over the duplicates and choose the lowest of the counts
names <- unique(dups$fullname)
for (i in 1:length(names)) {
  dups[which(dups$fullname == names[i]),1] <- min(dups[which(dups$fullname == names[i]),1]) # replace all counts for taxon i with the minimum value for that taxon
}
# remove duplicates from dups
dups <- dups[!duplicated(dups[,6]),]
# put everything back together, remove indicator column and resolved_name column
n_cts_clean <- rbind(not_dups,dups)
n_cts_clean <- data.frame(cbind(n_cts_clean$fullname,n_cts_clean$minimum))
colnames(n_cts_clean) <- c("name","ct")

# remove lycophytes (CCDB classifies them as 'pteridophytes' so they come with the dataset)
n_cts_clean <- n_cts_clean[which(!n_cts_clean$name == "Selaginella_moellendorffii" & !n_cts_clean$name == "Lycopodium_clavatum"),]

#### Pruning the tree
## prune tree to count data
# load tree
tr <- ftolr::ftol_ml_dated_tree
# prune
tr_pruned <- drop.tip(tr, tr$tip.label[-na.omit(match(n_cts_clean$name, tr$tip.label))])
tr_pruned$tip.label
# write the tree
write.tree(tr_pruned, file = "data_subsets/polypodiopsida/polypodiopsida_csomepruned.tree") # write the pruned tree

#### Pruning the count data
# use pruned tree names to remove species not in the tree. 
tree_cts <- n_cts_clean[n_cts_clean$name %in% tr_pruned$tip.label,] # prune counts to tree


#### format data for analysis
## make final data frames for exporting
## make final data frames for exporting
# change the data frame to a *named vector* to make final .tsv for non-hidden state analyses
tree.counts <- as.numeric(tree_cts[,2])
names(tree.counts) <- tree_cts[,1]
pterid.counts <- data.frame(tree.counts)

# expand data frame for hidden states
maximum.ct <- max(tree.counts) + 10 # max chromosomes + 10 more (Mayrose et al. 2010 recommendation)
second_column <- as.numeric(tree_cts[,2])
expand_data <- second_column+maximum.ct+1 # Expansion (leave room for '0' in hidden state )

# creating the double entry in each datum
counts.expanded <- paste0("(", second_column, " ", expand_data,")")
names(counts.expanded) <- tree_cts[,1]
pterid.counts.hid <- data.frame(counts.expanded)


#### Writing the datasets
# write simple data to .tsv
write.table(pterid.counts, 
            file = "data_subsets/polypodiopsida/polypodiopsida_ccdb_mincounts.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)

# write hidden state data to .tsv
write.table(pterid.counts.hid, 
            file = "data_subsets/polypodiopsida/polypodiopsida_ccdb_mincounts_expanded.tsv",
            sep = "\t",
            col.names = FALSE,
            quote = FALSE)
