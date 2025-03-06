library(purrr) # map
library(rstan) # monitor
library(data.table) # fread
library(dplyr) # select

# Define constants
cohort <- c("99-00", "01-02", "07-08", "09-10", "11-12", "13-14", "15-16")
seeds <- c("6734", "4880", "5916")
pop <- c("Total", 
  "6 - 11 years", 
  "12 - 19 years", 
  "20 - 65 years",
  "over 65 years", 
  "0 - 5 years")

# Initialize report
opt <- "MCSim/report_mcmc.log"
cat("", file = opt) 
sink(opt, append=TRUE)  
cat(paste0("Starting time: ", Sys.time()), "\n\n")
sink()

# Process cohort
for (j in seq(7)) {
  #
  sink(opt, append=TRUE)  
  cat("\n", paste(cohort[j], "\n", "====================================\n"))
  sink()
  #
  if(j == 1) { 
    k <- 4 
  } else if (j == 7) {
    k <- 6
  } else k <- 5
  #
  for (age in seq(k)) {
        if(j %in% c(1,2)){
    out <- paste0("outputs/gPYR_analytic_5met_", cohort[j], "_", pop[age], "_", seeds, ".out")
  } else if (j %in% c(3:4)){
    out <- paste0("outputs/gPYR_analytic_4met_", cohort[j], "_", pop[age], "_", seeds, ".out")
  } else if (j %in% c(5:7)){
    out <- paste0("outputs/gPYR_analytic_3met_", cohort[j], "_", pop[age], "_", seeds, ".out")
  }
  #
  sink(opt, append=TRUE)  
  cat("\n", paste(pop[age], "\n")) 
  pars_name <- c("Ve_CUriFPBA(1)", "Ve_CUri3PBA(1)", "Ve_CUriDCCAt(1)",
    "Ve_CUriDBCA(1)", "Ve_CUriDCCAc(1)", 
    "lnCYF_IngDose(1)", "lnPRM_IngDose(1)", "lnDLM_IngDose(1)",
    "V_lnCYF_IngDose(1)", "V_lnPRM_IngDose(1)", "V_lnDLM_IngDose(1)",
    "LnPrior", "LnData", "LnPosterior")
  if(j %in% c(3:4)) pars_name <- pars_name[-c(5)]
  if(j %in% c(5:7)) pars_name <- pars_name[-c(4 ,5, 8, 11)]
  data <- out |> map(fread) |> map(select, all_of(pars_name)) |> map(as.data.frame) 
  n_chains <- length(data)
  sample_number <- dim(data[[1]])[1]
  dim <- c(sample_number, n_chains, dim(data[[1]])[2])
  n_iter <- dim(data[[1]])[1]
  n_param <- dim(data[[1]])[2]
  x <- array(sample_number:(n_iter * n_chains * n_param), dim = dim)
  for (i in 1:n_chains) {
    x[, i, ] <- as.matrix(data[[i]][1:n_iter, ])
  }
  dimnames(x)[[3]] <- names(data[[1]])
  dim(x)
  mnt <- monitor(x[, , pars_name], digit = 6, print = T)
  sink()
  #
  save_directory <- "MCSim/saves"
  iters <- 1:dim(x)[1]
  rand.samp <- sample(iters, 1000)
  pyrs <- c(
    "lnCYF_IngDose(1)", "lnPRM_IngDose(1)", "lnDLM_IngDose(1)",
    "V_lnCYF_IngDose(1)", "V_lnPRM_IngDose(1)", "V_lnDLM_IngDose(1)"
  )
  if (j %in% seq(4)) { 
    lPsampsPBK <- x[rand.samp, , pyrs]
    dim(lPsampsPBK) <- c(3000,6)
    colnames(lPsampsPBK) <- pyrs
  }
  else if (j %in% c(5:7)){
    lPsampsPBK <- x[rand.samp, , pyrs[-c(3,6)]]
    dim(lPsampsPBK) <- c(3000,4)
    colnames(lPsampsPBK) <- pyrs[-c(3,6)]
  }
  file_name <- paste("mcsim-IPsamps_", pop[age], "_", 
    format(Sys.time(), "%Y-%m-%d"), ".RData", sep = "")
  save(lPsampsPBK, file = file.path(save_directory, cohort[j], file_name))
  }
  #
}

for (j in seq(7)) {
  
  if(j %in% c(1,2)){
     out <- paste0("outputs/gPYR_analytic_5met_", cohort[j], "_Total_5916.out")
  } else if (j %in% c(3:4)){
    out <- paste0("outputs/gPYR_analytic_4met_", cohort[j], "_Total_5916.out")
  } else if (j %in% c(5:7)){
    out <- paste0("outputs/gPYR_analytic_3met_", cohort[j], "_Total_5916.out")
  } 
  
  # Read last line from each chain
  X <- list()
  for(i in 1:3){
    setnames(
      X[[i]] <- fread(cmd=paste0("tail -1 '", out[i], "'")),
      names(fread(out[1], nrows = 0))
    )
  } 
  # Run final MCMC check
  MCMCfinal <- do.call(rbind, list(X[[3]], X[[2]], X[[1]]))
  MCMCfinal |> tail(1) |> write.table(file="MCMC.check.dat", row.names=F, sep="\t")  
  model <- "gPYR_analytic_ss.model"
  if(j %in% c(1:2)){
    input <- paste0("gPYR_analytic_5met_check_", cohort[j], "_Total", ".mcmc.in")
  } else if (j %in% c(3:4)){
    input <- paste0("gPYR_analytic_4met_check_", cohort[j], "_Total", ".mcmc.in")
  } else if (j %in% c(5:7)){
    input <- paste0("gPYR_analytic_3met_check_", cohort[j], "_Total", ".mcmc.in")
  }
  input 
  RMCSim::mcsim(model = model, input = input, dir = "MCSim")
  newname <- paste0("mcmc.check_", cohort[j])
  file.rename("mcmc.check.out", newname)
}


