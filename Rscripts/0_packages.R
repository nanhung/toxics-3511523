## Ensure required tools are available
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", dependencies = TRUE)
}

# Fresh install RMCSim from GitHub
if (requireNamespace("RMCSim", quietly = TRUE)) {
  remove.packages("RMCSim")
  remotes::install_github("nanhung/RMCSim")
  RMCSim::install_mcsim(version = "6.1.0")
}

# Optional: Compile model code (uncomment if needed)
# model_a <- "gPYR_analytic_ss.model"
# RMCSim::makemcsim(model = model_a, mxstep = 5000, dir = "MCSim")
# model_p <- "gPYR_pbk.model"
# RMCSim::makemcsim(model = model_p, mxstep = 5000, dir = "MCSim")

# List of required packages
pkgs <- c(
  "ggplot2", "scales", "dplyr", "survey", "magrittr", "foreach",
  "doParallel", "readxl", "foreign", "purrr", "rstan", "data.table",
  "ggpubr", "forcats", "tidyr", "tibble", "reshape2", "EnvStats",
  "cowplot", "ggdendro", "deSolve"
)

# Install missing packages efficiently
missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, dependencies = TRUE, quiet = TRUE)
}

