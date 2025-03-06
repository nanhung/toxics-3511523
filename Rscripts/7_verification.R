# Reproduce TK graph ----------------------------------------------------------
library(RMCSim)

# Load data
model <- "gPYR_pbk.model"

# 161.179 (uM) * 434 (ug/uM) / 70 (kg) / 1000 (ug/mg)
# 1 mg CYF =  0.27 mg FPBA
file.copy("MCSim/mcsim.gPYR_pbk.model.exe", "mcsim.gPYR_pbk.model.exe")
mcsim(model = model, input = "gPYR_pbk_oDLM001.in", dir = "MCSim")
oDLM_sim <- read.delim("sim.out", skip = 1)
mcsim(model = model, input = "gPYR_pbk_oPRM01.in", dir = "MCSim")
oPRM_sim <- read.delim("sim.out", skip = 1)
mcsim(model = model, input = "gPYR_pbk_oCPM01.in", dir = "MCSim")
oCPM_sim <- read.delim("sim.out", skip = 1)
mcsim(model = model, input = "gPYR_pbk_oCYF003.in", dir = "MCSim")
oCYF_sim <- read.delim("sim.out", skip = 1)
mcsim(model = model, input = "gPYR_pbk_iCYF160.in", dir = "MCSim")
iCYF_sim <- read.delim("sim.out", skip = 1)
mcsim(model = model, input = "gPYR_pbk_dCPM31.in", dir = "MCSim")
dCPM_sim <- read.delim("sim.out", skip = 1)

#
dat <- read.csv("data/TKData.csv")
comp <- "DLM"
expo <- "Oral"
meta <- "3PBA"
oDLM_3PBA <- subset(dat, Compound == comp & Dose == 0.01 & Exposure == expo
                & Metabolite == meta)
meta <- "DBCA"
oDLM_DBCA <- subset(dat, Compound == comp & Dose == 0.01 & Exposure == expo 
                & Metabolite == meta)
comp <- "PRM"
expo <- "Oral"
meta <- "cis-DCCA"
oPRM_DCCAc <- subset(dat, Compound == comp & Dose == 0.1 & Exposure == expo
                & Metabolite == meta)
meta <- "3PBA"
oPRM_3PBA <- subset(dat, Compound == comp & Dose == 0.1 & Exposure == expo
                & Metabolite == meta)
meta <- "trans-DCCA"
oPRM_DCCAt <- subset(dat, Compound == comp & Dose == 0.1& Exposure == expo
                & Metabolite == meta)
comp <- "CPM"
expo <- "Oral"
meta <- "cis-DCCA"
oCPM_DCCAc <- subset(dat, Compound == comp & Dose == 0.1 & Exposure == expo
                & Metabolite == meta)
meta <- "3PBA"
oCPM_3PBA <- subset(dat, Compound == comp & Dose == 0.1 & Exposure == expo
                & Metabolite == meta)
meta <- "trans-DCCA"
oCPM_DCCAt <- subset(dat, Compound == comp & Dose == 0.1 & Exposure == expo
                & Metabolite == meta)
comp <- "CYF"
expo <- "Oral"
meta <- "cis-DCCA"
oCYF_DCCAc <- subset(dat, Compound == comp & Dose == 0.03 & Exposure == expo
                & Metabolite == meta)
meta <- "FPBA"
oCYF_FPBA <- subset(dat, Compound == comp & Dose == 0.03 & Exposure == expo
                & Metabolite == meta)
meta <- "trans-DCCA"
oCYF_DCCAt <- subset(dat, Compound == comp & Dose == 0.03 & Exposure == expo
                & Metabolite == meta)
comp <- "CYF"
expo <- "Inhalation"
meta <- "cis-DCCA"
iCYF_DCCAc <- subset(dat, Compound == comp & Dose == 0.000160 & Exposure == expo
                & Metabolite == meta)
meta <- "FPBA"
iCYF_FPBA <- subset(dat, Compound == comp & Dose == 0.000160 & Exposure == expo
                & Metabolite == meta)
meta <- "trans-DCCA"
iCYF_DCCAt <- subset(dat, Compound == comp & Dose == 0.000160 & Exposure == expo
                & Metabolite == meta)
comp <- "CPM"
expo <- "Dermal"
meta <- "cis-DCCA"
dCPM_DCCAc <- subset(dat, Compound == comp & Dose == 31 & Exposure == expo
                & Metabolite == meta)
meta <- "3PBA"
dCPM_3PBA <- subset(dat, Compound == comp & Dose == 31 & Exposure == expo
                & Metabolite == meta)
meta <- "trans-DCCA"
dCPM_DCCAt <- subset(dat, Compound == comp & Dose == 31 & Exposure == expo
                & Metabolite == meta)

#
library(dplyr)
library(tibble)

verification_ratio <- function(pred = oDLM_sim, obs = oDLM_3PBA, variable = "CUri3PBA") {
  pred |> filter(Time %in% obs$Time) |> select(all_of(variable)) |> 
    add_column(observed = obs$value) |> 
    `colnames<-`(c("predicted", "observed")) |>
    mutate(ratio = predicted/observed) }
verification_ratio()

x <- do.call(rbind, 
  list(verification_ratio(), 
    verification_ratio(pred = oDLM_sim, obs = oDLM_DBCA, variable = "CUriDBCA"), 
    verification_ratio(pred = oDLM_sim, obs = oDLM_3PBA, variable = "CUri3PBA"), 
    verification_ratio(pred = oPRM_sim, obs = oPRM_DCCAc, variable = "CUriDCCAc"),
    verification_ratio(pred = oPRM_sim, obs = oPRM_DCCAt, variable = "CUriDCCAt"),
    verification_ratio(pred = oPRM_sim, obs = oPRM_3PBA, variable = "CUri3PBA"),
    verification_ratio(pred = oCPM_sim, obs = oCPM_3PBA, variable = "CUri3PBA"),
    verification_ratio(pred = oCPM_sim, obs = oCPM_DCCAc, variable = "CUriDCCAc"),
    verification_ratio(pred = oCPM_sim, obs = oCPM_DCCAt, variable = "CUriDCCAt"),
    verification_ratio(pred = oCYF_sim, obs = oCYF_DCCAc, variable = "RAMetDCCAc"),
    verification_ratio(pred = oCYF_sim, obs = oCYF_DCCAt, variable = "RAMetDCCAt"),
    verification_ratio(pred = oCYF_sim, obs = oCYF_FPBA, variable = "RAMetFPBA"),
    verification_ratio(pred = iCYF_sim, obs = iCYF_FPBA, variable = "RAMetFPBA"),
    verification_ratio(pred = iCYF_sim, obs = iCYF_DCCAc, variable = "RAMetDCCAc"),
    verification_ratio(pred = iCYF_sim, obs = iCYF_DCCAt, variable = "RAMetDCCAt"),
    verification_ratio(pred = dCPM_sim, obs = dCPM_DCCAc, variable = "RAMetDCCAc"),
    verification_ratio(pred = dCPM_sim, obs = dCPM_DCCAt, variable = "RAMetDCCAt"),
    verification_ratio(pred = dCPM_sim, obs = dCPM_3PBA, variable = "RAMet3PBA")
    )
)
x |> dim()

x |> mutate(accept = ifelse(ratio > 10, 0, ifelse(ratio < 0.1, 0, 1))) |> 
  summarise(rate = sum(accept)/n())
x |> mutate(accept = ifelse(ratio > 3, 0, ifelse(ratio < 0.333, 0, 1))) |> 
  summarise(rate = sum(accept)/n())
x |> mutate(accept = ifelse(ratio > 2, 0, ifelse(ratio < 0.5, 0, 1))) |> 
  summarise(rate = sum(accept)/n())

library(reshape2)
df_sim <- do.call(rbind, 
  list(
    melt(oDLM_sim, "Time") |> add_column(Compound = "DLM", Exposure = "Oral"),
    melt(oPRM_sim, "Time") |> add_column(Compound = "PRM", Exposure = "Oral"),
    melt(oCPM_sim, "Time") |> add_column(Compound = "CPM", Exposure = "Oral"),
    melt(oCYF_sim, "Time") |> add_column(Compound = "CYF", Exposure = "Oral"),
    melt(iCYF_sim, "Time") |> add_column(Compound = "CYF", Exposure = "Inhalation"),
    melt(dCPM_sim, "Time") |> add_column(Compound = "CPM", Exposure = "Dermal")))
dim(df_sim)

library(ggplot2)
custom_theme <- theme(strip.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(colour = "black"),
                      plot.title.position='plot',
                      legend.text = element_text(size=11),
                      strip.text = element_text(size=11),
                      axis.title = element_text(size=10),
                      axis.text = element_text(size=10))

p1 <- melt(oDLM_sim, "Time") |> 
  mutate(Metabolite = ifelse(variable=="CUri3PBA", "3PBA", "DBCA")) |>
  ggplot() + 
  geom_line(aes(x=Time, y = value)) +
  facet_wrap(~Metabolite) +
  geom_point(data = oDLM_3PBA, aes(x= Time, y=value)) +
  geom_point(data = oDLM_DBCA, aes(x= Time, y=value)) +
  scale_y_log10() +
  ggtitle("Oral deltamethrin (Sam et al., 2012)") +
  xlab("Time (h)") +
  ylab("Urinary concentration (mg/l)") +
  theme_bw() + custom_theme
p2 <- melt(oPRM_sim, "Time") |> 
  mutate(Metabolite = ifelse(variable=="CUri3PBA", "3PBA", 
      ifelse(variable=="CUriDCCAc", "cis-DCCA", "trans-DCCA"))) |>
  ggplot() + 
  geom_line(aes(x=Time, y = value)) +
  facet_wrap(~Metabolite) +
  geom_point(data = oPRM_3PBA, aes(x=Time, y=value)) +
  geom_point(data = oPRM_DCCAc, aes(x=Time, y=value)) +
  geom_point(data = oPRM_DCCAt, aes(x=Time, y=value)) +
  scale_y_log10() +
  ggtitle("Oral permethrin (Ratelle et al., 2015)") +
  xlab("Time (h)") +
  ylab("Urinary concentration (mg/l)") +
  theme_bw() + custom_theme
p3 <- melt(oCPM_sim, "Time") |> 
  mutate(Metabolite = ifelse(variable=="CUri3PBA", "3PBA", 
      ifelse(variable=="CUriDCCAc", "cis-DCCA", "trans-DCCA"))) |>
  ggplot() + 
  geom_line(aes(x=Time, y = value)) +
  facet_wrap(~Metabolite) +
  geom_point(data = oCPM_3PBA, aes(x=Time, y=value)) +
  geom_point(data = oCPM_DCCAc, aes(x=Time, y=value)) +
  geom_point(data = oCPM_DCCAt, aes(x=Time, y=value)) +
  scale_y_log10() +
  ggtitle("Oral cypermethrin (Ratelle et al., 2015)") +
  xlab("Time (h)") +
  ylab("Urinary concentration (mg/l)") +
  theme_bw() + custom_theme
p4 <- melt(oCYF_sim, "Time") |> 
  mutate(Metabolite = ifelse(variable=="RAMetFPBA", "FPBA", 
      ifelse(variable=="RAMetDCCAc", "cis-DCCA", "trans-DCCA"))) |>
  ggplot() + 
  geom_line(aes(x=Time, y = value)) +
  facet_wrap(~Metabolite) +
  geom_point(data = oCYF_FPBA, aes(x=Time, y=value)) +
  geom_point(data = oCYF_DCCAc, aes(x=Time, y=value)) +
  geom_point(data = oCYF_DCCAt, aes(x=Time, y=value)) +
  scale_y_log10() +
  ggtitle("Oral cyfluthrin (Leng et al., 1997)") +
  xlab("Time (h)") +
  ylab("Ecxretion rate (ug/h)") +
  theme_bw() + custom_theme
p5 <- melt(iCYF_sim, "Time") |> 
  mutate(Metabolite = ifelse(variable=="RAMetFPBA", "FPBA", 
      ifelse(variable=="RAMetDCCAc", "cis-DCCA", "trans-DCCA"))) |>
  ggplot() + 
  geom_line(aes(x=Time, y = value)) +
  facet_wrap(~Metabolite) +
  geom_point(data = iCYF_FPBA, aes(x=Time, y=value)) +
  geom_point(data = iCYF_DCCAc, aes(x=Time, y=value)) +
  geom_point(data = iCYF_DCCAt, aes(x=Time, y=value)) +
  scale_y_log10() +
  ggtitle("Inhalation cyfluthrin (Leng et al., 1997)") +
  xlab("Time (h)") +
  ylab("Ecxretion rate (ug/h)") +
  theme_bw() + custom_theme
p5 <- melt(dCPM_sim, "Time") |> 
  mutate(Metabolite = ifelse(variable=="RAMet3PBA", "3PBA", 
      ifelse(variable=="RAMetDCCAc", "cis-DCCA", "trans-DCCA"))) |>
  ggplot() + 
  geom_line(aes(x=Time, y = value)) +
  facet_wrap(~Metabolite) +
  geom_point(data = dCPM_3PBA, aes(x=Time, y=value)) +
  geom_point(data = dCPM_DCCAc, aes(x=Time, y=value)) +
  geom_point(data = dCPM_DCCAt, aes(x=Time, y=value)) +
  scale_y_log10() +
  ggtitle("Dermal cypermethrin (Woollen et al., 1992)") +
  xlab("Time (h)") +
  ylab("Ecxretion rate (ug/h)") +
  theme_bw() + custom_theme

library(cowplot)
plot_grid(p1, p2, p3, p4, p5, ncol = 1)
ggsave("verification.pdf", width = 10, height=12)

# housekeeping
file.remove("sim.out")
file.remove("mcsim.gPYR_pbk.model.exe")
