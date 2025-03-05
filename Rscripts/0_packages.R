# Install RMCSim package from personal repo
if (!require(remotes)) install.packages("remotes")
if (require(RMCSim)) remove.packages("RMCSim")
remotes::install_github("nanhung/RMCSim")
RMCSim::install_mcsim(version = "6.1.0")

# Compile model code (optional, the model program is compiled in MCSim folder)
#model_a <- "gPYR_analytic_ss.model"
#RMCSim::makemcsim(model = model_a, mxstep = 5000, dir = "MCSim")
#model_p <- "gPYR_pbk.model"
#RMCSim::makemcsim(model = model_p, mxstep = 5000, dir = "MCSim")

# Install other R packages
pkgs <- c("ggplot2", "scales", "dplyr", "survey",  "magrittr", "foreach",
  "doParallel", "readxl", "foreign", "purrr", "rstan", "data.table", "ggpubr",
  "forcats", "tidyr", "tibble", "reshape2", "EnvStats", "cowplot", "ggdendro", "deSolve")

for (i in seq(length(pkgs))){
  pac <- pkgs[i]
  if (!require(pac)) install.packages(pac)
}
