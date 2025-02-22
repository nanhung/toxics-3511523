# Install & Load package -------------------------------------------------------
# Use Linux (or WSL)

#if (require(RMCSim)) remove.packages("RMCSim")
#remotes::install_github("nanhung/RMCSim")

library(RMCSim)
#install_mcsim(version = "6.1.0")
#mcsim_version()

# Compile model code -----------------------------------------------------------
model <- "gPYR_analytic_ss.model"
#model <- "gPYR.model"
#makemcsim(model = model, mxstep = 5000, dir = "MCSim")

# Choose input file ------------------------------------------------------------
#input <- "gPYR_analytic.mcmc.in"

# MCMC -------------------------------------------------------------------------
library(foreach)
library(doParallel)

system("rm *.out")
#system("rm -rf outputs")
if (!dir.exists("outputs")) dir.create("outputs")

current_files <- list.files()
detectCores()
cores <- 3    
cl <- makeCluster(cores)
registerDoParallel(cl)

strt <- Sys.time()
out <- foreach(i = 1:cores) %dopar% {
  set.seed(i + 1)
  RMCSim::mcsim(model = model, input = input, dir = "MCSim", parallel = T, check=F)
}
print(Sys.time() - strt)

new_files <- setdiff(list.files(), current_files)
to_remove <- new_files[grep(".kernel|.in", new_files)]
file.remove(to_remove)
out_files <- setdiff(list.files(), current_files)

for (i in 1:cores) file.copy(out_files[i], paste0("outputs/", out_files[i]))
#file.remove(out_files)

