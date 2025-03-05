library(RMCSim)
library(dplyr)

# Load model
file.copy("MCSim/mcsim.gPYR_pbk.model.exe", "mcsim.gPYR_pbk.model.exe")
model <- "gPYR_PBPK.model"

# AUC calculation
auc <- function(dat){
  time_interval <- 0.05
  n_obs <- length(dat)
  sum(time_interval*(dat[-1] + dat[-n_obs])/2)
}

# Deltamethrin
mcsim(model = model, input = "gPYR_6y_DLM.mtc.in", dir = "MCSim")
mtc_sim_dlm <- read.delim("gPYR_6y_DLM.mtc.out")
CPLS <- grep("CPLS", colnames(mtc_sim_dlm))
CBRN <- grep("CBRN", colnames(mtc_sim_dlm))
sim_dlm_cpls <- mtc_sim_dlm[ , CPLS]
sim_dlm_cbrn <- mtc_sim_dlm[ , CBRN]

dlm_df <- data.frame(
  compound = "Deltamethrin",
  metrics = c(rep("Oral dose", 100),
    rep("Plasma Cmax", 100), rep("Brain Cmax", 100), 
    rep("Plasma AUC", 100), rep("Brain AUC", 100)),
  moe = c(
    1.79 / mtc_sim_dlm[, grep("IngDose_DLM", colnames(mtc_sim_dlm))],
    118.2 / apply(sim_dlm_cpls, 1, max), 
    34 / apply(sim_dlm_cbrn, 1, max),
    1174.1 / apply(sim_dlm_cpls, 1, auc),
    475 / apply(sim_dlm_cbrn, 1, auc))
)
dlm_df |> group_by(metrics) |>
  summarise(P50 = median(moe), P05=quantile(moe, .05))

# Permethrin
mcsim(model = model, input = "gPYR_1y_PRM.mtc.in", dir = "MCSim")
mtc_sim_prm <- read.delim("gPYR_1y_PRM.mtc.out")
CPM_CPLS <- grep("CPLS_5", colnames(mtc_sim_prm))
CPM_CBRN <- grep("CBRN_5", colnames(mtc_sim_prm))
sim_cpm_cpls <- mtc_sim_prm[ , CPM_CPLS]
sim_cpm_cbrn <- mtc_sim_prm[ , CPM_CBRN]
TPM_CPLS <- grep("CPLS_6", colnames(mtc_sim_prm))
TPM_CBRN <- grep("CBRN_6", colnames(mtc_sim_prm))
sim_tpm_cpls <- mtc_sim_prm[ , TPM_CPLS]
sim_tpm_cbrn <- mtc_sim_prm[ , TPM_CBRN]

cpm_df <- data.frame(
  compound = "cis-Permethrin",
  metrics = c(
    rep("Oral dose", 100),
    rep("Plasma Cmax", 100), rep("Brain Cmax", 100), 
    rep("Plasma AUC", 100), rep("Brain AUC", 100)),
  moe = c( 
    39.17 / mtc_sim_prm[, grep("IngDose_PRM", colnames(mtc_sim_prm))],
    1605.5 / apply(sim_cpm_cpls, 1, max), 
    459.8 / apply(sim_cpm_cbrn, 1, max),
    15822.6 / apply(sim_cpm_cpls, 1, auc),
    6400 / apply(sim_cpm_cbrn, 1, auc))
)
cpm_df |> group_by(metrics) |> 
  summarise(P50 = median(moe), P05=quantile(moe, .05))

tpm_df <- data.frame(
  compound = "trans-Permethrin",
  metrics = c(
    rep("Oral dose", 100),
    rep("Plasma Cmax", 100), rep("Brain Cmax", 100), 
    rep("Plasma AUC", 100), rep("Brain AUC", 100)),
  moe = c(
    39.17 / mtc_sim_prm[, grep("IngDose_PRM", colnames(mtc_sim_prm))],
    503 / apply(sim_cpm_cpls, 1, max), 
    143.8 / apply(sim_cpm_cbrn, 1, max),
    4926.1 / apply(sim_tpm_cpls, 1, auc),
    1996.4 / apply(sim_tpm_cbrn, 1, auc))
)
tpm_df |> group_by(metrics) |> summarise(P50 = median(moe), P05=quantile(moe, .05))


# Cyfluthrin
mcsim(model = model, input = "gPYR_1y_CYF.mtc.in", dir = "MCSim")
mtc_sim_cyf <- read.delim("gPYR_1y_CYF.mtc.out")
CPLS_CCF <- grep("CPLS_1", colnames(mtc_sim_cyf))
CBRN_CCF <- grep("CBRN_1", colnames(mtc_sim_cyf))
sim_cpls_ccf <- mtc_sim_cyf[ , CPLS_CCF]
sim_cbrn_ccf <- mtc_sim_cyf[ , CBRN_CCF]
CPLS_TCF <- grep("CPLS_2", colnames(mtc_sim_cyf))
CBRN_TCF <- grep("CBRN_2", colnames(mtc_sim_cyf))
sim_cpls_tcf <- mtc_sim_cyf[ , CPLS_TCF]
sim_cbrn_tcf <- mtc_sim_cyf[ , CBRN_TCF]

ccf_df <- data.frame(
  compound = "cis-Cyfluthrin",
  metrics = c(
    rep("Oral dose", 100),
    rep("Plasma Cmax", 100), rep("Brain Cmax", 100), 
    rep("Plasma AUC", 100), rep("Brain AUC", 100)),
  moe = c(
     1.22 / mtc_sim_cyf[, grep("IngDose_CYF", colnames(mtc_sim_cyf))],
     32.2 / apply(sim_cpls_ccf, 1, max), 
     9.26 / apply(sim_cbrn_ccf, 1, max),
     319.8 / apply(sim_cpls_ccf, 1, auc),
     129.4 / apply(sim_cbrn_ccf, 1, auc))
)
ccf_df |> group_by(metrics) |> 
  summarise(P50 = median(moe), P05=quantile(moe, .05))

tcf_df <- data.frame(
  compound = "trans-Cyfluthrin",
  metrics = c(
    rep("Oral dose", 100),
    rep("Plasma Cmax", 100), rep("Brain Cmax", 100), 
    rep("Plasma AUC", 100), rep("Brain AUC", 100)),
  moe = c(
     1.22 / mtc_sim_cyf[, grep("IngDose_CYF", colnames(mtc_sim_cyf))],
    48.3 / apply(sim_cpls_tcf, 1, max), 
    13.9 / apply(sim_cbrn_tcf, 1, max),
    479.8 / apply(sim_cpls_tcf, 1, auc),
    194.1 / apply(sim_cbrn_tcf, 1, auc))
)
tcf_df |> group_by(metrics) |> 
  summarise(P50 = median(moe), P05=quantile(moe, .05))

# Cypermethrin
mcsim(model = model, input = "gPYR_1y_CPM.mtc.in", dir = "MCSim")
mtc_sim_cpm <- read.delim("gPYR_1y_CPM.mtc.out")
CPLS_CCP <- grep("CPLS_3", colnames(mtc_sim_cpm))
CBRN_CCP <- grep("CBRN_3", colnames(mtc_sim_cpm))
sim_cpls_ccp <- mtc_sim_cpm[ , CPLS_CCP]
sim_cbrn_ccp <- mtc_sim_cpm[ , CBRN_CCP]
CPLS_TCP <- grep("CPLS_4", colnames(mtc_sim_cpm))
CBRN_TCP <- grep("CBRN_4", colnames(mtc_sim_cpm))
sim_cpls_tcp <- mtc_sim_cpm[ , CPLS_TCP]
sim_cbrn_tcp <- mtc_sim_cpm[ , CBRN_TCP]

ccp_df <- data.frame(
  compound = "cis-Cypermethrin",
  metrics = c(
    rep("Oral dose", 100),
    rep("Plasma Cmax", 100), rep("Brain Cmax", 100), 
    rep("Plasma AUC", 100), rep("Brain AUC", 100)),
  moe = c(
     5.48 / mtc_sim_cpm[, grep("IngDose_CPM", colnames(mtc_sim_cpm))],
     152.2 / apply(sim_cpls_ccp, 1, max), 
     1510.8 / apply(sim_cbrn_ccp, 1, max),
     43.7 / apply(sim_cpls_ccp, 1, auc),
     611.2 / apply(sim_cbrn_ccp, 1, auc))
)
ccp_df |> group_by(metrics) |> 
  summarise(P50 = median(moe), P05=quantile(moe, .05))

tcp_df <- data.frame(
  compound = "trans-Cypermethrin",
  metrics = c(
    rep("Oral dose", 100),
    rep("Plasma Cmax", 100), rep("Brain Cmax", 100), 
    rep("Plasma AUC", 100), rep("Brain AUC", 100)),
  moe = c( 
    5.48 / mtc_sim_cpm[, grep("IngDose_CPM", colnames(mtc_sim_cpm))],
    210.4 / apply(sim_cpls_tcp, 1, max), 
    60.5 / apply(sim_cbrn_tcp, 1, max),
    2087.9 / apply(sim_cpls_tcp, 1, auc),
    845.7 / apply(sim_cbrn_tcp, 1, auc))
)
tcp_df |> group_by(metrics) |> 
  summarise(P50 = median(moe), P05=quantile(moe, .05))

# Plot
library(ggplot2)
library(forcats)
library(scales)

metrics_level <- c("Oral dose", "Brain Cmax", "Brain AUC", "Plasma Cmax",
                   "Plasma AUC")
compound_level <- c("Deltamethrin",
                    "cis-Cyfluthrin", "trans-Cyfluthrin",
                    "cis-Permethrin", "trans-Permethrin",
                    "cis-Cypermethrin", "trans-Cypermethrin")
tot_df <- rbind(dlm_df, tpm_df, cpm_df, ccf_df, tcf_df, ccp_df, tcp_df) |>
  mutate(across(metrics, as.factor)) |>
  mutate(metrics = fct_relevel(metrics, metrics_level)) |>
  mutate(across(compound, as.factor)) |>
  mutate(compound = fct_relevel(compound, compound_level))

png("MOE.png", width = 2500, height = 1500, res = 300)
tot_df |> ggplot(aes(x = compound, y = moe, fill = compound)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_violin() +
  ylab("MOE/MOIE") +
  facet_grid(~metrics, switch = "x") +
  scale_fill_brewer() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    strip.text = element_text(size = 11),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11),
  )
dev.off()

file.remove("mcsim.gPYR_pbk.model.exe")
