library(dplyr)
library(ggplot2)
library(scales)

to_read <- c(
  "MCSim/saves/MCMC.check_9900.out",
  "MCSim/saves/MCMC.check_0102.out",
  "MCSim/saves/MCMC.check_0708.out",
  "MCSim/saves/MCMC.check_0910.out",
  "MCSim/saves/MCMC.check_1112.out",
  "MCSim/saves/MCMC.check_1314.out",
  "MCSim/saves/MCMC.check_1516.out")
cohort <- c("99-00", "01-02", "07-08", "09-10", "11-12", "13-14", "15-16")

for (i in seq(7)){
  x <- read.delim(to_read[i])
  x$cohort <- cohort[i]
  if (i == 1) X <- x else X <- rbind(X, x)
}

#
df <- X |> 
  mutate(Var = ifelse(Output_Var == "CUri3PBA", "3PBA", 
    ifelse(Output_Var == "CUriFPBA", "FPBA", 
      ifelse(Output_Var == "CUriDBCA", "DBCA", 
        ifelse(Output_Var == "CUriDCCAc", "cis-DCCA", 
          ifelse(Output_Var == "CUriDCCAt", "trans-DCCA", "<LOD")))))) 

df$Var <- factor(df$Var, level = c("FPBA", "3PBA", "trans-DCCA", 
    "cis-DCCA", "DBCA", "<LOD"))
df$cohort <- factor(df$cohort, level = cohort)

set_theme <- theme(
  legend.position  = c(0.5, 0.15),
  axis.text        = element_text(color = "black", size = 13),
  axis.ticks.y     = element_line(color = "black"),
  axis.ticks.x     = element_line(color = "black"),
  axis.line.x      = element_line(color = "black"),
  axis.line.y      = element_line(color = "black"),
  legend.key       = element_blank(),
  strip.background = element_blank(),
  axis.title       = element_text(size = 13),
  strip.text       = element_text(size = 13),
  legend.text      = element_text(size=12),
  legend.title     = element_blank(),  
  panel.background = element_blank()
)

p1 <- ggplot(data = df, aes(x=Data, y=Prediction, color = Var, shape = Var)) + 
  geom_point(alpha=0.5) + 
  scale_x_log10(#limits=c(0.1,1000),
    breaks = trans_breaks("log10", function(x) 10^x, n = 4),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(#limits=c(0.1,1000),
    breaks = trans_breaks("log10", function(x) 10^x, n = 4),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_colour_viridis_d(end = 0.95 ) +
  geom_abline(slope = 1) + 
  xlab("NHANES biomonitoring data (ug/L)") +
  ylab("PBK prediction (ug/L)") +
  facet_wrap(~cohort) +
  #geom_smooth(method='lm', se=F) +
  geom_abline(slope = 1, intercept = 0.477, linetype = "dashed", color = "grey") +
  geom_abline(slope = 1, intercept = -0.477, linetype = "dashed", color = "grey") +
  #ggtitle("Comparison of PBPK model prediction and NHANES biomonitoring data") +
  theme_bw() +
  set_theme +
  xlab("NHANES biomonitoring data (ug/L)") +
  ylab("PBPK prediction (ug/L)")

png("fig3_goodness.png", width = 3200, height = 2700, res = 350)
p1
dev.off()
