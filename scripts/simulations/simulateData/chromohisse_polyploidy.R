setwd("/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/2.Methods/simulation/src")
source("simulate.R", chdir = TRUE)

############
# settings #
############

# states 1:num_chromo are chromosome numbers in hidden state A
# states num_chromo + 1:num_chromo are chromosome numbers in hidden state B

time       <- 50
num_chromo <- 100 # must be even and >10
num_hidden <- 2
num_cats   <- num_hidden * num_chromo

######################
# specify parameters #
######################

# adjusted values for simulations
gamma_a <- 0.00081
delta_a <- 0.00765
rho_a   <- 0.00315
eta_a   <- 0.00585

gamma_b <- 0.00162
delta_b <- 0.0153
rho_b   <- 0.0063
eta_b   <- 0.0117

q_ab <- 0.0117
q_ba <- 0.0117

clado_no_change_a   <- 0.61650
clado_fission_a     <- 0.00068
clado_fusion_a      <- 0.00495
clado_polyploidy_a  <- 0.03105
clado_demipoly_a    <- 0.00585

clado_no_change_b   <- 0.61650
clado_fission_b     <- 0.00136
clado_fusion_b      <- 0.0099
clado_polyploidy_b  <- 0.0621
clado_demipoly_b    <- 0.0117

clado_ab <- 0.0315
clado_ba <- 0.0315
mu <- 0.549

absolute_clado_rates  <- c(clado_no_change_a, clado_no_change_b, clado_fission_a, clado_fission_b, clado_fusion_a, clado_fusion_b, clado_polyploidy_a, clado_polyploidy_b, clado_demipoly_a, clado_demipoly_b, clado_ab)
total_speciation_rate <- sum(absolute_clado_rates)
relative_clado_rates  <- absolute_clado_rates / total_speciation_rate

###################################
# spontaneous state-change events #
###################################

# spontaneous state changes
H = matrix(0, num_cats, num_cats)

# hidden state A
for(i in 1:num_chromo) { # i is the state
  
  # i is the index of state i in hidden state A
  i_hid = i + num_chromo # i_hid is the index of state i in hidden state B
  
  if (i == 1) { # n=1 (minimum chromosome number): fission, polyploidy, demipolyploidy
    H[i, i + 1] <- gamma_a + rho_a + 0.5 * eta_a 
    H[i_hid, i_hid + 1] <- gamma_b + rho_b + 0.5 * eta_b 
    
  } else if (i == 2) { # n=2: fission, fusion, polyploidy (+2), and even demipolyploidy (+1)
    H[i, i - 1] <- delta_a
    H[i_hid, i_hid - 1] <- delta_b
    H[i, i + 1] <- gamma_a + eta_a
    H[i_hid, i_hid + 1] <- gamma_b + eta_b
    H[i, i + i] <- rho_a
    H[i_hid, i_hid + i] <- rho_b
    
  } else if (i == 3) { # n=3: fission, fusion, polyploidy (+i), and odd demipolyploidy (+i*1.5)
    H[i, i - 1] <- delta_a
    H[i_hid, i_hid - 1] <- delta_b
    H[i, i + 1] <- gamma_a + 0.5 * eta_a
    H[i_hid, i_hid + 1] <- gamma_b + 0.5 * eta_b
    H[i, i + i] <- rho_a
    H[i_hid, i_hid + i] <- rho_b
    H[i, i + 2] <- 0.5 * eta_a
    H[i_hid, i_hid + 2] <- 0.5 * eta_b
    
  } else if (3 < i & i + i <= num_chromo) { # 3 < i < k/2: fission, fusion, polyploidy, & demipolyploidy
    H[i, i - 1] <- delta_a
    H[i_hid, i_hid - 1] <- delta_b
    H[i, i + 1] <- gamma_a
    H[i_hid, i_hid + 1] <- gamma_b
    H[i, i + i] <- rho_a
    H[i_hid, i_hid + i] <- rho_b
    
    if (i %% 2 == 0) { # even demipolyploidy
      H[i, i * 1.5] <- eta_a
      H[i_hid, i * 1.5 + num_chromo] <- eta_b
      
    } else if (i %% 2 != 0) { # odd demipolyploidy
      H[i, i * 1.5 - 0.5] <- 0.5 * eta_a
      H[i_hid, i * 1.5 - 0.5 + num_chromo] <- 0.5 * eta_b
      H[i, i * 1.5 + 0.5] <- 0.5 * eta_a
      H[i_hid, i * 1.5 + 0.5 + num_chromo] <- 0.5 * eta_b
    }
    
  } else if (num_chromo < i + i & i * 1.5 <= num_chromo) { # k/2 < i < k/1.5: # fission, fusion, & demipolyploidy
    H[i, i - 1] <- delta_a
    H[i_hid, i_hid - 1] <- delta_b
    H[i, i + 1] <- gamma_a
    H[i_hid, i_hid + 1] <- gamma_b
    
    if (i %% 2 == 0) { # even demipolyploidy
      H[i, i * 1.5] <- eta_a
      H[i_hid, i * 1.5 + num_chromo] <- eta_b
      
    } else if (i %% 2 != 0) { # odd demipolyploidy
      H[i, i * 1.5 - 0.5] <- 0.5 * eta_a
      H[i_hid, i * 1.5 - 0.5 + num_chromo] <- 0.5 * eta_b
      H[i, i * 1.5 + 0.5] <- 0.5 * eta_a
      H[i_hid, i * 1.5 + 0.5 + num_chromo] <- 0.5 * eta_b
    }
    
  } else if (i * 1.5 - 0.5 == num_chromo) { # boundary case when odd-descending, but not odd-ascending, demipolyploidy is allowed: # fission, fusion, & reducing demipolyploidy
    H[i, i - 1] <- delta_a
    H[i_hid, i_hid - 1] <- delta_b
    H[i, i + 1] <- gamma_a
    H[i_hid, i_hid + 1] <- gamma_b
    H[i, i * 1.5 - 0.5] <- 0.5 * eta_a
    H[i_hid, i * 1.5 - 0.5 + num_chromo] <- 0.5 * eta_b
    
  } else if (i < num_chromo) { # fission & fusion only
    H[i, i - 1] <- delta_a
    H[i_hid, i_hid - 1] <- delta_b
    H[i, i + 1] <- gamma_a
    H[i_hid, i_hid + 1] <- gamma_b
    
  } else if (i == num_chromo) { # n=k (max chromosome number): fusion only
    H[i, i - 1] <- delta_a
    H[i_hid, i_hid - 1] <- delta_b
  }

  H[i, i_hid] <- q_ab
  H[i_hid, i] <- q_ba
}


diag(H) = -rowSums(H)

##################################################
# diversification-associated state-change events #
##################################################

# empty speciation-rate functions
lambda <- numeric(num_cats)

# empty speciation-associated state-change array
Omega = array(0, dim = c(num_cats, num_cats, num_cats))

for(i in 1:num_chromo) {
  
  i_hid = i + num_chromo 
  
  if (i == 1) { # fissions, demipolyploidy, polyploidy 
  
    # Hidden State A
    lambda[i] <- clado_no_change_a + clado_fission_a + clado_polyploidy_a + 0.5*clado_demipoly_a + clado_ab
    # no change
    Omega[i,i,i]          <- 1.0 * clado_no_change_a / lambda[i]
    # increase by 1
    Omega[i,i,i + 1]      <- 0.5 * (clado_fission_a + clado_polyploidy_a + (0.5*clado_demipoly_a)) / lambda[i]
    Omega[i,i + 1,i]      <- 0.5 * (clado_fission_a + clado_polyploidy_a + (0.5*clado_demipoly_a)) / lambda[i]
    # hidden state change
    Omega[i,i,i_hid]      <- 0.5 * clado_ab / lambda[i]
    Omega[i,i_hid,i]      <- 0.5 * clado_ab / lambda[i]
    
    # Hidden State B
    lambda[i_hid] <- clado_no_change_b + clado_fission_b + clado_polyploidy_b + 0.5*clado_demipoly_b + clado_ba
    # no change
    Omega[i_hid,i_hid,i_hid]        <- 1.0 * clado_no_change_b / lambda[i_hid]
    # increase by 1
    Omega[i_hid,i_hid,i_hid + 1]    <- 0.5 * (clado_fission_b + clado_polyploidy_b + (0.5*clado_demipoly_b)) / lambda[i_hid]
    Omega[i_hid,i_hid + 1,i_hid]    <- 0.5 * (clado_fission_b + clado_polyploidy_b + (0.5*clado_demipoly_b)) / lambda[i_hid]
    # hidden state change
    Omega[i_hid,i_hid,i]            <- 0.5 * clado_ba / lambda[i_hid]
    Omega[i_hid,i,i_hid]            <- 0.5 * clado_ba / lambda[i_hid]
    
  } else if (i == 2) { # fissions, fusions, demipolyploidy, polyploidy
  
    # Hidden State A
    lambda[i] <- clado_no_change_a + clado_fusion_a + clado_fission_a + clado_polyploidy_a + clado_demipoly_a + clado_ab
    # no change
    Omega[i,i,i]          <- 1.0 * clado_no_change_a / lambda[i]
    # increase by 1
    Omega[i,i,i + 1]      <- 0.5 * (clado_fission_a + clado_demipoly_a) / lambda[i]
    Omega[i,i + 1,i]      <- 0.5 * (clado_fission_a + clado_demipoly_a) / lambda[i]
    # decrease by 1
    Omega[i,i,i - 1]      <- 0.5 * clado_fusion_a  / lambda[i]
    Omega[i,i - 1,i]      <- 0.5 * clado_fusion_a  / lambda[i]
    # polyploidy
    Omega[i,i,i*2]        <- 0.5 * clado_polyploidy_a  / lambda[i]
    Omega[i,i*2,i]        <- 0.5 * clado_polyploidy_a  / lambda[i]
    # hidden state change
    Omega[i,i,i_hid]      <- 0.5 * clado_ab / lambda[i]
    Omega[i,i_hid,i]      <- 0.5 * clado_ab / lambda[i]
    
    # Hidden State B
    lambda[i_hid] <- clado_no_change_b + clado_fusion_b + clado_fission_b + clado_polyploidy_b + clado_demipoly_b + clado_ba
    # no change
    Omega[i_hid,i_hid,i_hid]        <- 1.0 * clado_no_change_b / lambda[i_hid]
    # increase by 1
    Omega[i_hid,i_hid,i_hid + 1]    <- 0.5 * (clado_fission_b + clado_demipoly_b) / lambda[i_hid]
    Omega[i_hid,i_hid + 1,i_hid]    <- 0.5 * (clado_fission_b + clado_demipoly_b) / lambda[i_hid]
    # decrease by 1
    Omega[i_hid,i_hid,i_hid - 1]    <- 0.5 * clado_fusion_b  / lambda[i_hid]
    Omega[i_hid,i_hid - 1,i_hid]    <- 0.5 * clado_fusion_b  / lambda[i_hid]
    # polyploidy
    Omega[i_hid,i_hid,i*2+num_chromo]      <- 0.5 * clado_polyploidy_b  / lambda[i_hid]
    Omega[i_hid,i*2+num_chromo,i_hid]      <- 0.5 * clado_polyploidy_b  / lambda[i_hid]
    # hidden state change
    Omega[i_hid,i_hid,i]            <- 0.5 * clado_ba / lambda[i_hid]
    Omega[i_hid,i,i_hid]            <- 0.5 * clado_ba / lambda[i_hid]
    
  } else if (i == 3) { # fissions, fusions, demipolyploidy, polyploidy
  
    # Hidden State A
    lambda[i] <- clado_no_change_a + clado_fusion_a + clado_fission_a + clado_polyploidy_a + clado_demipoly_a + clado_ab
    # no change
    Omega[i,i,i]          <- 1.0 * clado_no_change_a / lambda[i]
    # increase by 1
    Omega[i,i,i + 1]      <- 0.5 * (clado_fission_a + (0.5*clado_demipoly_a)) / lambda[i]
    Omega[i,i + 1,i]      <- 0.5 * (clado_fission_a + (0.5*clado_demipoly_a)) / lambda[i]
    # decrease by 1
    Omega[i,i,i - 1]      <- 0.5 * clado_fusion_a  / lambda[i]
    Omega[i,i - 1,i]      <- 0.5 * clado_fusion_a  / lambda[i]
    # polyploidy
    Omega[i,i,i*2]        <- 0.5 * clado_polyploidy_a  / lambda[i]
    Omega[i,i*2,i]        <- 0.5 * clado_polyploidy_a  / lambda[i]
    # demipolyploidy
    # up
    Omega[i,i,i+2]        <- 0.5 * clado_demipoly_a  / lambda[i]
    Omega[i,i+2,i]        <- 0.5 * clado_demipoly_a  / lambda[i]
    # hidden state change
    Omega[i,i,i_hid]      <- 0.5 * clado_ab / lambda[i]
    Omega[i,i_hid,i]      <- 0.5 * clado_ab / lambda[i]
    
    # Hidden State B
    lambda[i_hid] <- clado_no_change_b + clado_fusion_b + clado_fission_b + clado_polyploidy_b + clado_demipoly_b + clado_ba
    # no change
    Omega[i_hid,i_hid,i_hid]        <- 1.0 * clado_no_change_b / lambda[i_hid]
    # increase by 1
    Omega[i_hid,i_hid,i_hid + 1]    <- 0.5 * (clado_fission_b + (0.5*clado_demipoly_b)) / lambda[i_hid]
    Omega[i_hid,i_hid + 1,i_hid]    <- 0.5 * (clado_fission_b + (0.5*clado_demipoly_b)) / lambda[i_hid]
    # decrease by 1
    Omega[i_hid,i_hid,i_hid - 1]    <- 0.5 * clado_fusion_b  / lambda[i_hid]
    Omega[i_hid,i_hid - 1,i_hid]    <- 0.5 * clado_fusion_b  / lambda[i_hid]
    # polyploidy
    Omega[i_hid,i_hid,i*2+num_chromo]      <- 0.5 * clado_polyploidy_b  / lambda[i_hid]
    Omega[i_hid,i*2+num_chromo,i_hid]      <- 0.5 * clado_polyploidy_b  / lambda[i_hid]
    # demipolyploidy
    # up
    Omega[i_hid,i_hid,i_hid+2]      <- 0.5 * clado_demipoly_b  / lambda[i_hid]
    Omega[i_hid,i_hid+2,i_hid]      <- 0.5 * clado_demipoly_b  / lambda[i_hid]
    # hidden state change
    Omega[i_hid,i_hid,i]            <- 0.5 * clado_ba / lambda[i_hid]
    Omega[i_hid,i,i_hid]            <- 0.5 * clado_ba / lambda[i_hid]
    
  } else if (3 < i && i + i <= num_chromo) { # fissions, fusions, demipolyploidy, polyploidy
  
    # Hidden State A
    lambda[i] <- clado_no_change_a + clado_fusion_a + clado_fission_a + clado_polyploidy_a + clado_demipoly_a + clado_ab
    # no change
    Omega[i,i,i]          <- 1.0 * clado_no_change_a / lambda[i]
    # increase by 1
    Omega[i,i,i + 1]      <- 0.5 * clado_fission_a / lambda[i]
    Omega[i,i + 1,i]      <- 0.5 * clado_fission_a / lambda[i]
    # decrease by 1
    Omega[i,i,i - 1]      <- 0.5 * clado_fusion_a  / lambda[i]
    Omega[i,i - 1,i]      <- 0.5 * clado_fusion_a  / lambda[i]
    # polyploidy
    Omega[i,i,i*2]        <- 0.5 * clado_polyploidy_a  / lambda[i]
    Omega[i,i*2,i]        <- 0.5 * clado_polyploidy_a  / lambda[i]
    # demipolyploidy
    if (i %% 2 == 0) { # even chromo counts
      Omega[i,i,i*1.5]    <- 0.5 * clado_demipoly_a  / lambda[i]
      Omega[i,i*1.5,i]    <- 0.5 * clado_demipoly_a  / lambda[i]
    } else { # odd chromo counts
      Omega[i,i,i*1.5+0.5] <- 0.25 * clado_demipoly_a  / lambda[i]
      Omega[i,i*1.5+0.5,i] <- 0.25 * clado_demipoly_a  / lambda[i]
      Omega[i,i,i*1.5-0.5] <- 0.25 * clado_demipoly_a  / lambda[i]
      Omega[i,i*1.5-0.5,i] <- 0.25 * clado_demipoly_a  / lambda[i]
    }
    # hidden state change
    Omega[i,i,i_hid]      <- 0.5 * clado_ab / lambda[i]
    Omega[i,i_hid,i]      <- 0.5 * clado_ab / lambda[i]
    
    # Hidden State B
    lambda[i_hid] <- clado_no_change_b + clado_fusion_b + clado_fission_b + clado_polyploidy_b + clado_demipoly_b + clado_ba
    # no change
    Omega[i_hid,i_hid,i_hid]          <- 1.0 * clado_no_change_b / lambda[i_hid]
    # increase by 1
    Omega[i_hid,i_hid,i_hid + 1]      <- 0.5 * clado_fission_b / lambda[i_hid]
    Omega[i_hid,i_hid + 1,i_hid]      <- 0.5 * clado_fission_b / lambda[i_hid]
    # decrease by 1
    Omega[i_hid,i_hid,i_hid - 1]      <- 0.5 * clado_fusion_b  / lambda[i_hid]
    Omega[i_hid,i_hid - 1,i_hid]      <- 0.5 * clado_fusion_b  / lambda[i_hid]
    # polyploidy
    Omega[i_hid,i_hid,i*2+num_chromo]        <- 0.5 * clado_polyploidy_b  / lambda[i_hid]
    Omega[i_hid,i*2+num_chromo,i_hid]        <- 0.5 * clado_polyploidy_b  / lambda[i_hid]
    # demipolyploidy
    if (i %% 2 == 0) { # even chromo counts
      Omega[i_hid,i_hid,i*1.5+num_chromo]    <- 0.5 * clado_demipoly_b  / lambda[i_hid]
      Omega[i_hid,i*1.5+num_chromo,i_hid]    <- 0.5 * clado_demipoly_b  / lambda[i_hid]
    } else { # odd chromo counts
      Omega[i_hid,i_hid,i*1.5+0.5+num_chromo]  <- 0.25 * clado_demipoly_b  / lambda[i_hid]
      Omega[i_hid,i*1.5+0.5+num_chromo,i_hid]  <- 0.25 * clado_demipoly_b  / lambda[i_hid]
      Omega[i_hid,i_hid,i*1.5-0.5+num_chromo]  <- 0.25 * clado_demipoly_b  / lambda[i_hid]
      Omega[i_hid,i*1.5-0.5+num_chromo,i_hid]  <- 0.25 * clado_demipoly_b  / lambda[i_hid]
    }
    # hidden state change
    Omega[i_hid,i_hid,i]              <- 0.5 * clado_ba / lambda[i_hid]
    Omega[i_hid,i,i_hid]              <- 0.5 * clado_ba / lambda[i_hid]
    
  } else if (num_chromo < i + i &&  i * 1.5 <= num_chromo) { # fissions, fusions, demipolyploidy
  
    # Hidden State A
    lambda[i] <- clado_no_change_a + clado_fusion_a + clado_fission_a + clado_demipoly_a + clado_ab
    # no change
    Omega[i,i,i]          <- 1.0 * clado_no_change_a / lambda[i]
    # increase by 1
    Omega[i,i,i + 1]      <- 0.5 * clado_fission_a / lambda[i]
    Omega[i,i + 1,i]      <- 0.5 * clado_fission_a / lambda[i]
    # decrease by 1
    Omega[i,i,i - 1]      <- 0.5 * clado_fusion_a  / lambda[i]
    Omega[i,i - 1,i]      <- 0.5 * clado_fusion_a  / lambda[i]
    # demipolyploidy
    if (i %% 2 == 0) { # even chromo counts
      Omega[i,i,i*1.5]    <- 0.5 * clado_demipoly_a  / lambda[i]
      Omega[i,i*1.5,i]    <- 0.5 * clado_demipoly_a  / lambda[i]
    } else { # odd chromo counts
      Omega[i,i,i*1.5+0.5] <- 0.25 * clado_demipoly_a  / lambda[i]
      Omega[i,i*1.5+0.5,i] <- 0.25 * clado_demipoly_a  / lambda[i]
      Omega[i,i,i*1.5-0.5] <- 0.25 * clado_demipoly_a  / lambda[i]
      Omega[i,i*1.5-0.5,i] <- 0.25 * clado_demipoly_a  / lambda[i]
    }
    # hidden state change
    Omega[i,i,i_hid]       <- 0.5 * clado_ab / lambda[i]
    Omega[i,i_hid,i]       <- 0.5 * clado_ab / lambda[i]
    
    # Hidden State B
    lambda[i_hid] <- clado_no_change_b + clado_fusion_b + clado_fission_b + clado_demipoly_b + clado_ba
    # no change
    Omega[i_hid,i_hid,i_hid]          <- 1.0 * clado_no_change_b / lambda[i_hid]
    # increase by 1
    Omega[i_hid,i_hid,i_hid + 1]      <- 0.5 * clado_fission_b / lambda[i_hid]
    Omega[i_hid,i_hid + 1,i_hid]      <- 0.5 * clado_fission_b / lambda[i_hid]
    # decrease by 1
    Omega[i_hid,i_hid,i_hid - 1]      <- 0.5 * clado_fusion_b  / lambda[i_hid]
    Omega[i_hid,i_hid - 1,i_hid]      <- 0.5 * clado_fusion_b  / lambda[i_hid]
    # demipolyploidy
    if (i %% 2 == 0) { # even chromo counts
      Omega[i_hid,i_hid,i*1.5+num_chromo]    <- 0.5 * clado_demipoly_b  / lambda[i_hid]
      Omega[i_hid,i*1.5+num_chromo,i_hid]    <- 0.5 * clado_demipoly_b  / lambda[i_hid]
    } else { # odd chromo counts
      Omega[i_hid,i_hid,i*1.5+0.5+num_chromo] <- 0.25 * clado_demipoly_b  / lambda[i_hid]
      Omega[i_hid,i*1.5+0.5+num_chromo,i_hid] <- 0.25 * clado_demipoly_b  / lambda[i_hid]
      Omega[i_hid,i_hid,i*1.5-0.5+num_chromo] <- 0.25 * clado_demipoly_b  / lambda[i_hid]
      Omega[i_hid,i*1.5-0.5+num_chromo,i_hid] <- 0.25 * clado_demipoly_b  / lambda[i_hid]
    }
    # hidden state change
    Omega[i_hid,i_hid,i]              <- 0.5 * clado_ba / lambda[i_hid]
    Omega[i_hid,i,i_hid]              <- 0.5 * clado_ba / lambda[i_hid]
    
  } else if (i * 1.5 - 0.5 == num_chromo) { # fissions, fusions, descending demipolyploidy
  
    # Hidden State A
    lambda[i] <- clado_no_change_a + clado_fusion_a + clado_fission_a + 0.5*clado_demipoly_a + clado_ab
    # no change
    Omega[i,i,i]          <- 1.0 * clado_no_change_a / lambda[i]
    # increase by 1
    Omega[i,i,i + 1]      <- 0.5 * clado_fission_a / lambda[i]
    Omega[i,i + 1,i]      <- 0.5 * clado_fission_a / lambda[i]
    # decrease by 1
    Omega[i,i,i - 1]      <- 0.5 * clado_fusion_a  / lambda[i]
    Omega[i,i - 1,i]      <- 0.5 * clado_fusion_a  / lambda[i]
    # demipolyploidy
    Omega[i,i,i*1.5-0.5]  <- 0.25 * clado_demipoly_a  / lambda[i]
    Omega[i,i*1.5-0.5,i]  <- 0.25 * clado_demipoly_a  / lambda[i]
    # hidden state change
    Omega[i,i,i_hid]      <- 0.5 * clado_ab / lambda[i]
    Omega[i,i_hid,i]      <- 0.5 * clado_ab / lambda[i]
    
    # Hidden State B
    lambda[i_hid] <- clado_no_change_b + clado_fusion_b + clado_fission_b + 0.5*clado_demipoly_b + clado_ba
    # no change
    Omega[i_hid,i_hid,i_hid]          <- 1.0 * clado_no_change_b / lambda[i_hid]
    # increase by 1
    Omega[i_hid,i_hid,i_hid + 1]      <- 0.5 * clado_fission_b / lambda[i_hid]
    Omega[i_hid,i_hid + 1,i_hid]      <- 0.5 * clado_fission_b / lambda[i_hid]
    # decrease by 1
    Omega[i_hid,i_hid,i_hid - 1]      <- 0.5 * clado_fusion_b  / lambda[i_hid]
    Omega[i_hid,i_hid - 1,i_hid]      <- 0.5 * clado_fusion_b  / lambda[i_hid]
    # demipolyploidy
    Omega[i_hid,i_hid,i*1.5-0.5+num_chromo]  <- 0.25 * clado_demipoly_b  / lambda[i_hid]
    Omega[i_hid,i*1.5-0.5+num_chromo,i_hid]  <- 0.25 * clado_demipoly_b  / lambda[i_hid]
    # hidden state change
    Omega[i_hid,i_hid,i]              <- 0.5 * clado_ba / lambda[i_hid]
    Omega[i_hid,i,i_hid]              <- 0.5 * clado_ba / lambda[i_hid]
    
  } else if (i < num_chromo) { # fissions, fusions
  
    # Hidden State A
    lambda[i] <- clado_no_change_a + clado_fusion_a + clado_fission_a + clado_ab
    # no change
    Omega[i,i,i]          <- 1.0 * clado_no_change_a / lambda[i]
    # increase by 1
    Omega[i,i,i + 1]      <- 0.5 * clado_fission_a / lambda[i]
    Omega[i,i + 1,i]      <- 0.5 * clado_fission_a / lambda[i]
    # decrease by 1
    Omega[i,i,i - 1]      <- 0.5 * clado_fusion_a  / lambda[i]
    Omega[i,i - 1,i]      <- 0.5 * clado_fusion_a  / lambda[i]
    # hidden state change
    Omega[i,i,i_hid]      <- 0.5 * clado_ab / lambda[i]
    Omega[i,i_hid,i]      <- 0.5 * clado_ab / lambda[i]
    
    # Hidden State B
    lambda[i_hid] <- clado_no_change_b + clado_fusion_b + clado_fission_b + clado_ba
    # no change
    Omega[i_hid,i_hid,i_hid]      <- 1.0 * clado_no_change_b / lambda[i_hid]
    # increase by 1
    Omega[i_hid,i_hid,i_hid + 1]  <- 0.5 * clado_fission_b / lambda[i_hid]
    Omega[i_hid,i_hid + 1,i_hid]  <- 0.5 * clado_fission_b / lambda[i_hid]
    # decrease by 1
    Omega[i_hid,i_hid,i_hid - 1]  <- 0.5 * clado_fusion_b  / lambda[i_hid]
    Omega[i_hid,i_hid - 1,i_hid]  <- 0.5 * clado_fusion_b  / lambda[i_hid]
    # hidden state change
    Omega[i_hid,i_hid,i]          <- 0.5 * clado_ba / lambda[i_hid]
    Omega[i_hid,i,i_hid]          <- 0.5 * clado_ba / lambda[i_hid]
    
  } else if (i == num_chromo) { # fusion only
  
    # Hidden State A
    lambda[i] <- clado_no_change_a + clado_fusion_a + clado_ab
    # no change
    Omega[i,i,i]          <- 1.0 * clado_no_change_a / lambda[i]
    # decrease by 1
    Omega[i,i,i - 1]      <- 0.5 * clado_fusion_a  / lambda[i]
    Omega[i,i - 1,i]      <- 0.5 * clado_fusion_a  / lambda[i]
    # hidden state change
    Omega[i,i,i_hid]      <- 0.5 * clado_ab / lambda[i]
    Omega[i,i_hid,i]      <- 0.5 * clado_ab / lambda[i]
    
    # Hidden State B
    lambda[i_hid] <- clado_no_change_b + clado_fusion_b + clado_ba
    # no change
    Omega[i_hid,i_hid,i_hid]      <- 1.0 * clado_no_change_b / lambda[i_hid]
    # decrease by 1
    Omega[i_hid,i_hid,i_hid - 1]  <- 0.5 * clado_fusion_b  / lambda[i_hid]
    Omega[i_hid,i_hid - 1,i_hid]  <- 0.5 * clado_fusion_b  / lambda[i_hid]
    # hidden state change
    Omega[i_hid,i_hid,i]          <- 0.5 * clado_ba / lambda[i_hid]
    Omega[i_hid,i,i_hid]          <- 0.5 * clado_ba / lambda[i_hid]
  }
}


#############################################
# other asynchronous diversification events #
#############################################

extinction <- rep(mu, num_cats)

############
# sampling #
############

rho <- 0.15

############
# simulate #
############

for (j in 1:10) {
  
  # init <- num_chromo / 3
  init <- 15
  
  min_percent <- 0.1 # fraction of samples in each hidden state
  target_min <- 200  # minimum number of tips
  target_max <- 600  # maximum number of tips
  
  sim <- simulateChromoHiSSEConditional(lambda, extinction, H, Omega, rho, time, init, 
                                        num_chromo, num_hidden, min_percent, 
                                        target_min, target_max, verbose = TRUE)
  
  plot(sim)
  #sim$data
  
  
  ### write data to nexus 
  data <- data.frame(name = names(sim$dat), a = sim$data)
  data$a[data$a>num_chromo] <- data$a[data$a>num_chromo] - num_chromo
  data$b <- data$a + num_chromo + 1 # add one to leave room for a 0 state in the hidden
  data$new <- paste0("(", data$a, " ",data$b,")")
  data <- data[,c(1,4)]
  
  datfilename <- paste0("sims/sim", j, ".tsv")
  write.table(data, 
              file=datfilename,
              quote = FALSE,
              row.names = FALSE,
              sep="\t",
              col.names = FALSE)
  ### write tree to file
  treefilename <- paste0("sims/sim", j, ".tre")
  ape::write.tree(sim, 
                  file = treefilename)
}
