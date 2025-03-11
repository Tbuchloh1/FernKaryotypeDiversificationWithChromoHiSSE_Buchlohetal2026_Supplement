library(RevGadgets)

workingdir <- "/Users/tbuchloh/Dropbox/2.Dissertation/Projects/1.KaryotypeEvol_Ferns/2.Methods/tests/sim_reduced"
# Set working directory
setwd(workingdir)

# create empty dataframe
df <- data.frame(matrix(nrow = 100, ncol = 22))
colnames(df) <- c("sim", "alpha","chi", 
                  "clado_demipoly_a", "clado_demipoly_b", 
                  "clado_fission_a", "clado_fission_b", 
                  "clado_fusion_a", "clado_fusion_b", 
                  "clado_no_change_a", "clado_no_change_b", 
                  "clado_polyploid_a", "clado_polyploid_b", 
                  "delta_a", "delta_b", 
                  "eta_a", "eta_b", 
                  "gamma_a", "gamma_b", 
                  "rho_a", "rho_b",
                  "turnover")
df$sim <- 1:100
sim_vals <- c(0, 1.03E-02, 2.42E-02, # sim# & alpha & chi
              5.82E-03, 5.82E-02, # clado_demipoly_a & _b
              6.71E-04, 3.27E-02, # clado_fission_a & _b
              2.95E-03, 2.51E-02, # clado_fussion_a & _b
              5.82E-01, 2.42E-02, # clado_no_change_a & _b
              3.09E-02, 2.24E-01, # clado_polyploid_a & _b
              5.37E-03, 2.86E+00, # delta_a & _b
              4.92E-03, 1.21E-01, # eta_a & _b
              6.71E-04, 3.36E-01, # gamma_a & _b
              2.77E-03, 2.19E-01, # rho_a & _b
              4.92E-01) # turnover

for (i in 1:10) {
  #load traces
  chisse_path <- paste0("output/chisse_chisse/chisse_chisse_",i,".log")
  trace_1 <- readTrace(paths = chisse_path, burnin = 2000)
  trace_1 <- trace_1[[1]]

  #calculate summary stats
  for (j in 2:22) {
    col <- colnames(df)[j]
    med <- median(trace_1[,which(colnames(trace_1) == col)])
    stdev <- sd(trace_1[,which(colnames(trace_1) == col)])
    df[i,j] <- (med - sim_vals[j]) / stdev
  }
}


# df_norm <- sweep(df, 2, sim_vals, "-")

colors <- c("grey", "darkgrey", # sim# & alpha & chi
            "pink", "red", # clado_demipoly_a & _b
            "orange", "darkorange", # clado_fission_a & _b
            "lightyellow", "yellow", # clado_fussion_a & _b
            "lightgreen", "green", # clado_no_change_a & _b
            "lightblue", "blue", # clado_polyploid_a & _b
            "lightyellow", "yellow", # delta_a & _b
            "pink", "red", # eta_a & _b
            "orange", "darkorange", # gamma_a & _b
            "lightblue", "blue", # rho_a & _b
            "purple") # turnover
par(mar = c(9, 4, 4, 2))
boxplot(df[1:10,2:22], las = 2, ylim = c(-13,2), main = "Deviance of Model Estimate from Simulation Value: 
        ((median of posterior distribution) - (simulated value)) / (st.dev of posterior distribution)", col = colors, ylab = "Normalized, Scaled Estimates (simulations 1-100)")
abline(h = 0, lty = 2)


#how to scale these deviances from the median...
  # divide by the standard deviation

write.csv(df, file = "sim_summary.csv")
  
  


