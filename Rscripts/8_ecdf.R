library(ggplot2)
library(dplyr)
library(reshape2)
library(forcats)
library(scales)
load("nhanes_data.RData")

nhanes_data |> head()
nhanes_data$AgeGroup <- as.character(nhanes_data$AgeGroup)
nhanes_data$AgeGroup |> class()


png("figa1_ecdf.png", width = 3000, height = 1800, res = 250)
nhanes_data |> select(year, AgeGroup, URX4FP, URXOPM, URXTCC, URXCB3, URXCCC) |> 
  `colnames<-`(c("year", "AgeGroup", "FPBA", "3PBA", "tran-DCCA", "DBCA", "cis-DCCA")) |>
  mutate(
    cohort = ifelse(year == 1999, "99\u201300", 
      ifelse(year == 2001, "01\u201302", 
        ifelse(year == 2007, "07\u201308", 
          ifelse(year == 2009, "09\u201310",
            ifelse(year == 2011, "11\u201312",
              ifelse(year == 2013, "13\u201314",
                ifelse(year == 2015, "15\u201316", NA)))))))) |>
  melt() |> #group_by(year, AgeGroup, variable) |>
  mutate(cohort=fct_relevel(cohort, c("99\u201300", 
        "01\u201302", "07\u201308", "09\u201310", "11\u201312", "13\u201314", "15\u201316"))) |>
  mutate(AgeGroup = replace(AgeGroup, AgeGroup == "Over 65", "> 65")) |>
  mutate(AgeGroup = replace(AgeGroup, AgeGroup == "0 - 5", "< 6")) |>
  mutate(
    AgeGroup = ifelse(AgeGroup == "6 - 11", "6 \u2013 11", 
      ifelse(AgeGroup == "12 - 19", "12 \u2013 19", 
        ifelse(AgeGroup == "20 - 65", "20 \u2013 65", AgeGroup)))) |>
  mutate(AgeGroup=fct_relevel(AgeGroup, c("< 6", "6 \u2013 11", "12 \u2013 19", "20 \u2013 65", "> 65"))) |> 
  ggplot(aes(value, color=AgeGroup)) + stat_ecdf() +
  facet_grid(variable~cohort) +
  theme_bw() +
  theme(legend.position = "top",
    plot.title = element_text(face='bold'),
    strip.background=element_blank(),
    axis.title = element_text(size=14),
    axis.text = element_text(size=14),
    legend.text = element_text(size=13),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    strip.text = element_text(size=14)) +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n = 4),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  xlab("Urinary concentration (ng/ml)")
dev.off()
