# Install RMCSim package from personal repo
if (!require(remotes)) install.packages("remotes")
if (require(RMCSim)) remove.packages("RMCSim")
remotes::install_github("nanhung/RMCSim")
RMCSim::install_mcsim(version = "6.1.0")

# Compile model code (optional, the model program is compiled in MCSim folder)
#model <- "gPYR_analytic_ss.model"
#RMCSim::makemcsim(model = model, mxstep = 5000, dir = "MCSim")

# Install other R packages
pkgs <- c("ggplot2", "dplyr", "survey",  "magrittr", "foreach", "doParallel", 
          "readxl", "foreign")

for (i in seq(length(pkgs))){
  pac <- pkgs[i]
  if (!require(pac)) install.packages(pac)
}
