library(dplyr)      # Data manipulation
library(ggplot2)    # Plotting
library(scales)     # Scale transformations

# Define file paths and cohorts
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

# Read and combine data
df <- X |> 
  mutate(Var = ifelse(Output_Var == "CUri3PBA", "3PBA", 
    ifelse(Output_Var == "CUriFPBA", "FPBA", 
      ifelse(Output_Var == "CUriDBCA", "DBCA", 
        ifelse(Output_Var == "CUriDCCAc", "cis-DCCA", 
          ifelse(Output_Var == "CUriDCCAt", "trans-DCCA", "<LOD")))))) 

df$Var <- factor(df$Var, level = c("FPBA", "3PBA", "trans-DCCA", 
    "cis-DCCA", "DBCA", "<LOD"))
df$cohort <- factor(df$cohort, level = cohort)

# Define custom theme
set_theme <- theme(
  legend.position  = c(0.5, 0.1),
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

# Create plot
p1 <- ggplot(data = df, aes(x=Data, y=Prediction, color = Var, shape = Var)) + 
  geom_point(alpha=0.5, size=1) + 
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n = 4),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n = 4),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_colour_viridis_d(end = 0.95 ) +
  geom_abline(slope = 1) + 
   facet_wrap(~cohort) +
  geom_abline(slope = 1, intercept = 0.477, linetype = "dashed", color = "grey") +
  geom_abline(slope = 1, intercept = -0.477, linetype = "dashed", color = "grey") +
  theme_bw() +
  xlab("NHANES biomonitoring data (µg/L)") +
  ylab("PBK prediction (µg/L)") +
  set_theme

# Save plot
ggsave("fig3_goodness.png", p1, width = 7, height = 6, dpi = 300, units = "in")

