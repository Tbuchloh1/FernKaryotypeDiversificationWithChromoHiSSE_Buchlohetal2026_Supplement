library(Rcpp)
library(ape)
library(phangorn)
library(dispRity)

# load simulator
sourceCpp("simulate.cpp")

# simulate wrapper
simulateChromoHiSSEConditional <- function(lambda, mu, H, Omega, rho, time, initial_state,
                                    num_chromo, num_hidden, min_percent, target_min, target_max, verbose = FALSE) {
    
    # make sure the tree has a root at the start time
    repeat {
        
        tree <- simulateChromoHiSSE(lambda, mu, H, Omega, rho, time, initial_state, num_chromo, num_hidden, min_percent, target_min, target_max, verbose)  
        if ( abs(max(branching.times(tree)) - time) < 1e-6 ) {
            break
        } else {
            if (verbose) cat("rejected because root did not survive\n")
        }
        
    }
    
    return(tree)
    
    
}


simulateChromoHiSSE <- function(lambda, mu, H, Omega, rho, time, initial_state,
                                num_chromo, num_hidden, min_percent, target_min, target_max, verbose) {
    
    # turn Omega into a vector
    num_states <- length(lambda)
    omega <- vector("list", num_states)
    for(i in 1:num_states) {
        omega[[i]] <- Omega[i,,]
    }
    
    # data for which states belong to which hidden states
    state_to_hidden_state_map <- rep(1:num_hidden, each = num_chromo)

    # 0-index the initial state
    init <- initial_state - 1
        
    repeat {

        # cat("simulating...\n")
        
        # simulate the data
        sim <- simulateChromoHiSSEInternal(lambda, mu, H, omega, rho, time, init)
        
        # check if we hit max number of lineages
        if(is.null(sim)) {
            cat("max number of lineages reached")
            next
        }
        
        # reformat
        lineages <- data.frame(anc = sim$anc, desc = sim$desc, start_time = sim$start_time, end_time = sim$end_time, state = 1 + sim$state, status = sim$status)
        
        # make sure we have extant tips
        num_extant  <- sum(lineages$end_time == 0 & lineages$status == "sampled")
        
        if (verbose) cat(num_extant, "\n")
        
        success <- TRUE
        if (num_extant < target_min) {
            success <- FALSE
            if (verbose) cat("rejected because not enough species\n")
        }
        
        if (num_extant > target_max) {
            success <- FALSE
            if (verbose) cat("rejected because too many species\n")
        }
        
        # check that hidden states are all represented in at least x% of tips
        if (success) {
            
            # get the hidden states for each tip
            hidden_states_per_tip <- state_to_hidden_state_map[lineages$state[lineages$status == "sampled"]]
            
            # make sure each hidden state is represented
            if (any((tabulate(hidden_states_per_tip, nbins = num_hidden) / num_extant) < min_percent)) {
                if (verbose) cat("rejected because not enough tips in each hidden state\n")
                success <- FALSE
            }
            
        }
        
        if (success) {
            break
        }

    }
    
    # create an index map for phylo object    
    old_tip_indexes <- lineages$desc[lineages$status %in% c("extinct", "sampled")]
    new_tip_indexes <- 1:length(old_tip_indexes)
    
    old_root_index <- 0
    new_root_index <- length(old_tip_indexes) + 1
    
    old_node_indexes <- unique(lineages$anc)
    old_node_indexes <- old_node_indexes[!old_node_indexes %in% c(0)]
    new_node_indexes <- new_root_index + seq_along(old_node_indexes)
    
    index_map <- data.frame(old_index = c(old_tip_indexes, old_root_index, old_node_indexes),
                            new_index = c(new_tip_indexes, new_root_index, new_node_indexes))
    index_map <- index_map[order(index_map$old_index),]
    
    # create the phylo object with all the lineages
    edge <- as.matrix(lineages[,1:2])
    edge[,1] <- index_map$new_index[match(edge[,1], index_map$old_index)]
    edge[,2] <- index_map$new_index[match(edge[,2], index_map$old_index)]
    edge.length <- (lineages$start_time - lineages$end_time)
    tip.label <- paste0(ifelse(lineages[match(old_tip_indexes, lineages$desc),]$end_time == 0, "e", "f"), new_tip_indexes)
    Nnode <- length(new_node_indexes) + 1
    
    phylo <- list(edge = edge,
                  edge.length = edge.length,
                  tip.label = tip.label,
                  Nnode = Nnode)
    class(phylo) <- "phylo"
    
    # drop unsampled lineages
    extinct_tips    <- lineages$desc[lineages$status == "extinct"]
    indexes_to_drop <- index_map$new_index[match(extinct_tips, index_map$old_index)]
    tips_to_drop    <- phylo$tip.label[indexes_to_drop]
    pruned_tree     <- drop.tip(phylo, tips_to_drop)
    
    # create data matrix
    sampled_tips <- pruned_tree$tip.label
    sample_new_tip_indexes <- match(sampled_tips, tip.label)
    sample_old_tip_indexes <- index_map$old_index[match(sample_new_tip_indexes, index_map$new_index)]
    sampled_states <- lineages$state[match(sample_old_tip_indexes, lineages$desc)]
    names(sampled_states) <- sampled_tips
    
    # combine
    pruned_tree$data <- sampled_states
    
    return(ladderize(pruned_tree))
    
    
}