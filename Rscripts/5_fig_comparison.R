# Packages
library(dplyr) # bind_rows
library(ggpubr) # ggscatter; need cmake
library(scales)
library(forcats) # fct_reveal
library(tidyr)   # unite
library(tibble) # add_column
library(reshape2) # melt
library(EnvStats) # GeoMean
library(cowplot)  # plot_grid
library(ggdendro)

# 
all_cohorts_data <- c("99-00", "01-02", 
                      "07-08", "09-10", 
                      "11-12", "13-14", "15-16")

pyr_P <-c("68359-37-5", # CYF CAS
          "52645-53-1", # PRM CAS
          "52918-63-5", # DLM CAS
          "52315-07-8") # CPM CAS

for (cohort in all_cohorts_data){
  #
  OnlyPparms3 <- list.files(paste0("./bayesmarker/saves/", cohort,"/"), pattern='lPsamps-gm')
  OnlyPparms3
  PBK_files <- list.files(paste0("./MCSim/saves/", cohort,"/"), pattern='mcsim')
  PBK_files
  #
  op_gp_no <- c(
    grep("20 - 65", x = OnlyPparms3), 
    grep("12 - 19", x = OnlyPparms3), 
    grep("6 - 11", x = OnlyPparms3), 
    grep("66 years", x = OnlyPparms3), 
    grep("0 - 5", x = OnlyPparms3))
  pk_gp_no <- c(
    grep("20 - 65", x = PBK_files), 
    grep("12 - 19", x = PBK_files), 
    grep("6 - 11", x = PBK_files),
    grep("over 65", x = PBK_files), 
    grep("0 - 5", x = PBK_files))
  op_gps <- OnlyPparms3[op_gp_no]
  op_gps
  pk_gps <- PBK_files[pk_gp_no]
  pk_gps
  #
  if(cohort %in% all_cohorts_data[c(2:6)]){
    no_group <- 4
  } else if (cohort %in% all_cohorts_data[1]) {
    no_group <- 3
  } else if (cohort %in% all_cohorts_data[7]) {
    no_group <- 5
  }
  #
  for(i in seq(no_group)){
    #
    randsamp <- sample(1:3000, 300)
    #
    OnlyPparms3_to_be_load <- paste0("bayesmarker/saves/", cohort, "/", op_gps[i])
    OnlyPparms3_to_be_load
    load(OnlyPparms3_to_be_load)
    lPsampsgm |> dim()
    lPsampsgm <- lPsampsgm[randsamp, ]
    #
    PBK_file_to_be_load <- paste0("MCSim/saves/", cohort, "/", pk_gps[i])
    PBK_file_to_be_load
    load(PBK_file_to_be_load)
    lPsampsPBK |> dim()
    lPsampsPBK <- lPsampsPBK[randsamp, ]
    #
    if(cohort %in% all_cohorts_data[c(2:6)]){
      if (i == 1){
        age <- "20 - 65"
      } else if( i == 2 ){
        age <- "12 - 19"
      } else if ( i == 3 ){
        age <- "6 - 11"
      } else if ( i == 4 ){
        age <- "> 65"
      }
    } else if (cohort %in% all_cohorts_data[1]){
      if (i == 1){
        age <- "20 - 65"
      } else if( i == 2 ){
        age <- "12 - 19"
      } else if (i == 3 ) {
        age <- "6 - 11" 
      }
    } else if (cohort %in% all_cohorts_data[7]){
      if (i == 1){
        age <- "20 - 65"
      } else if( i == 2 ){
        age <- "12 - 19"
      } else if ( i == 3 ){
        age <- "6 - 11"
      } else if ( i == 4 ){
        age <- "> 65"
      } else if ( i == 5) {
        age <- "< 6"
      }
    }
    #
    app1_cyf <- data.frame(c(exp(lPsampsgm[,which(colnames(lPsampsgm) == pyr_P[1])])),
                           "Cyfluthrin", "High-throughput", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app1_prm <- data.frame(c(exp(lPsampsgm[,which(colnames(lPsampsgm) == pyr_P[2])])),
                           "Permethrin", "High-throughput", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app1_dlm <- data.frame(c(exp(lPsampsgm[,which(colnames(lPsampsgm) == pyr_P[3])])),
                           "Deltamethrin", "High-throughput", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app1_cpm <- data.frame(c(exp(lPsampsgm[,which(colnames(lPsampsgm) == pyr_P[4])])),
                           "Cypermethrin", "High-throughput", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app2_cyf <- data.frame(c(exp(lPsampsPBK[,"lnCYF_IngDose(1)"])),
                           "Cyfluthrin", "PBK", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app2_prm <- data.frame(c(exp(lPsampsPBK[,"lnPRM_IngDose(1)"])),
                           "Permethrin", "PBK", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app2_cpm <- data.frame(c(exp(lPsampsPBK[,"lnPRM_IngDose(1)"])*0.1),
                           "Cypermethrin", "PBK", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    if(cohort %in% all_cohorts_data[c(1:4)]){
      app2_dlm <- data.frame(c(exp(lPsampsPBK[,"lnDLM_IngDose(1)"])),
                             "Deltamethrin", "PBK", cohort, age) |>
        `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    }
    #
    pbk_pop_cyf <- data.frame(
      c(lPsampsPBK[,"lnCYF_IngDose(1)"]),
      c(lPsampsPBK[,"V_lnCYF_IngDose(1)"]),
      "Cyfluthrin", cohort, age) |>
      `colnames<-`(c("mean", "variance", "Compound", "Cohort", "Age"))
    pbk_pop_prm <- data.frame(
      c(lPsampsPBK[,"lnPRM_IngDose(1)"]),
      c(lPsampsPBK[,"V_lnPRM_IngDose(1)"]),
      "Permethrin", cohort, age) |>
      `colnames<-`(c("mean", "variance", "Compound", "Cohort", "Age"))
    if(cohort %in% all_cohorts_data[c(1:4)]){
      pbk_pop_dlm <- data.frame(
        c(lPsampsPBK[,"lnDLM_IngDose(1)"]),
        c(lPsampsPBK[,"V_lnDLM_IngDose(1)"]),
        "Deltamethrin", cohort, age) |>
      `colnames<-`(c("mean", "variance", "Compound", "Cohort", "Age"))
    }
    #
    if ( cohort=="99-00" & i==1 ){
      X <- bind_rows(app1_cyf, app2_cyf, 
                     app1_dlm, app2_dlm,
                     app1_prm, app2_prm,
                     app1_cpm, app2_cpm)
      X_pbk <- bind_rows(pbk_pop_cyf, pbk_pop_prm, pbk_pop_dlm)
    } else {
      if (cohort %in% all_cohorts_data[c(1:4)]){
        X <- bind_rows(X,
          app1_cyf, app2_cyf, 
          app1_dlm, app2_dlm,
          app1_prm, app2_prm,
          app1_cpm, app2_cpm)
        X_pbk <- bind_rows(X_pbk, pbk_pop_cyf, pbk_pop_prm, pbk_pop_dlm)
      } else {
        X <- bind_rows(X,
          app1_cyf, app2_cyf, 
          app1_dlm, 
          app1_prm, app2_prm,
          app1_cpm, app2_cpm)
        X_pbk <- bind_rows(X_pbk, pbk_pop_cyf, pbk_pop_prm)
      }   
    }
  }
}

X$Cohort <- factor(X$Cohort, level = all_cohorts_data)
X$Age <- factor(X$Age, level = c("< 6", "6 - 11", "12 - 19", "20 - 65", "> 65"))
X$Approach <- factor(X$Approach, level = c("PBK", "High-throughput"))

X.PBK <- X |> filter(Approach == 'PBK') |> group_by(Compound, Age, Cohort) |>
  summarise(PBK_median = median(Prediction),
            PBK_upper = quantile(Prediction, 0.975),
            PBK_lower = quantile(Prediction, 0.025))
X.Bayesmarker <- X |> filter(Approach == 'High-throughput') |> 
  group_by(Compound, Age, Cohort) |>
  summarise(Bayesmarker_median = median(Prediction), 
            Bayesmarker_upper = quantile(Prediction, 0.975),
            Bayesmarker_lower = quantile(Prediction, 0.025))
XX <- full_join(X.Bayesmarker, X.PBK, by = c('Compound', 'Age', 'Cohort')) 

XX |> mutate(ratio = Bayesmarker_median/PBK_median) |> 
  mutate(consistency = ifelse(ratio > 10, 0, ifelse(ratio < 0.1, 0, 1))) |>
  na.omit() |>
  group_by(Compound) |> 
  summarise(Located10fold = sum(consistency)/n())

png("fig4_comparison.png", width = 1800, height = 1200, res = 300)
XX |>
  ggscatter(x="PBK_median", y="Bayesmarker_median", 
    color = "Age", shape = "Compound", size = 2) +
  #geom_linerange(aes(ymin=Bayesmarker_lower, ymax=Bayesmarker_upper, color=Age)) +
  #geom_linerange(aes(xmin=PBK_lower, xmax=PBK_upper, color=Age)) +
  scale_x_log10(
    limit = c(1e-12, 1e-3),
    breaks = trans_breaks("log10", function(x) 10^x, n = 4),
    labels = trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_log10(
    limit = c(1e-12, 1e-3),
    breaks = trans_breaks("log10", function(x) 10^x, n = 4),
    labels = trans_format("log10", scales::math_format(10^.x))) + 
  geom_abline(slope=1) + 
  geom_abline(slope=1, intercept = 1, linetype = 2) + 
  geom_abline(slope=1, intercept = -1, linetype = 2) +
  geom_abline(slope=1, intercept = 0.477, linetype = 3) + 
  geom_abline(slope=1, intercept = -0.477, linetype = 3) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title=element_text(size=14),
    axis.text=element_text(size=14)) + 
  xlab("PBK model") + ylab("High-throughput model")
dev.off()


# Variability 
# use 3PBA due to the detection rate > 50% 
load("nhanes_data.RData")
nhanes_data_opm <- nhanes_data |> group_by(year, AgeGroup) |> 
  summarise(
    opm50 = quantile(URXOPM, c(0.5), na.rm=T),
    opm95 = quantile(URXOPM, c(0.95), na.rm=T),
    cv = sd(URXOPM, na.rm=T)/mean(URXOPM, na.rm=T),
    variance = var(log(URXOPM), na.rm=T)) |> 
  mutate(ratio = opm95/opm50) |> select(-c(opm50, opm95))

cohort <- c(rep("99-00", 3), 
      rep(c("01-02", "07-08", "09-10", "11-12", "13-14", "15-16"), each=4), 
      "15-16")
Age <- rep(c("6 - 11", "12 - 19", "20 - 65", "> 65"), 7)
nhanes_data_opm$Cohort <- cohort
nhanes_data_opm$Age <- c(Age[-4], "< 6")
nhanes_data_opm$Compound <- "3PBA"
#nhanes_data[, -c(1:2)] |> print(n=28)
#nhanes_data[, -c(1:2)] |> select(ratio) |> range()

x_pbk_var <- X_pbk |>
  group_by(Cohort, Compound, Age) |>
  summarise(n=n(), variance = mean(variance)) |> 
  select(-n)
x_pbk_ratio <- X_pbk |> 
  slice(rep(1:n(), each = 100)) |> # variability - 100 individuals
  mutate(lnConc = rnorm(n(), mean, variance^0.5)) |> 
  mutate(Conc = exp(lnConc)) |> 
  group_by(Cohort, Compound, Age) |> 
  summarise(n=n(), P50 = median(Conc), P95 = quantile(Conc, 0.95)) |> 
  mutate(ratio = P95/P50) |>
  select(-c(n, P50, P95))
x_diff <- x_pbk_ratio |> filter(Compound == "Permethrin") |>
  full_join(nhanes_data_opm[,c(6,8,7,5)], by = c('Age', 'Cohort')) |> 
  mutate(diff = ratio.y/ratio.x) 
x_diff$Age <- factor(x_diff$Age, 
  level = c("< 6", "6 - 11", "12 - 19", "20 - 65", "> 65"))
x_diff$Cohort <- factor(x_diff$Cohort, 
  level = c("99-00", "01-02", "07-08", "09-10", "11-12", "13-14", "15-16"))
x_pbk_join <- full_join(x_pbk_var, x_pbk_ratio,
  by = c("Compound", "Cohort", "Age")) 
x <- x_pbk_join |> filter(Compound == "Permethrin") |>
  rbind(nhanes_data_opm[,c(6,8,7,4,5)])

p1 <- x_diff |> 
  ggplot(aes(x=Cohort, y=Age, size=diff)) + geom_point() +
  scale_size_continuous(breaks = c(1,1.5,2))+
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title.x=element_text(size=14),
    axis.title.y=element_blank(),
    axis.text = element_text(size=13),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
  )
p2 <- x |> unite(Cohort_Age, c("Cohort", "Age"), remove=F) |>
  mutate(across(Age, as.factor)) |>
  mutate(Age=fct_relevel(Age, levels(X$Age))) |> 
  ggplot(aes(x=ratio, y=Compound, group=Cohort_Age, color=Age)) + 
  geom_line(color=1, linetype=2, linewidth=0.1) +
  geom_point() + 
  xlab("P95/P50") +
  theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.title.x=element_text(size=14),
    axis.title.y=element_blank(),
    axis.text = element_text(size=13),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"))

png("fig5_variability.png", width = 1800, height = 1200, res = 200)
plot_grid(p2, p1, nrow = 2, rel_heights = c(1/3,2/3),  align="v")
dev.off()

# Trend 
x_pbk_prm <- X |> 
  filter(Approach == 'PBK' & Compound == "Permethrin") |> 
  group_by(Approach, Age, Cohort) |>
  summarise(med = median(Prediction),
            upr = quantile(Prediction, 0.975),
            lor = quantile(Prediction, 0.025))
x_bayesmarker_prm <- X |> 
  filter(Approach == 'High-throughput' & Compound == "Permethrin") |> 
  group_by(Approach, Age, Cohort) |>
  summarise(med = median(Prediction), 
            upr = quantile(Prediction, 0.975),
            lor = quantile(Prediction, 0.025))

p2 <- rbind(x_bayesmarker_prm, x_pbk_prm) |> 
    ggplot(aes(x = Cohort, y = med, group=Approach, color = Approach)) +
    geom_point(size = 3) +
    geom_line() +
    #geom_linerange(aes(ymin = lor, ymax=upr)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
      labels = trans_format("log10", scales::math_format(10^.x))) + 
    facet_wrap(~Age, scale = "free_x", ncol=5) + 
    scale_colour_discrete(labels=c('PBK', 'High-throughput')) +
    ylab("Permetrhin intake rate (mg/kg BW/day)") +
    theme_bw() + 
    theme(
      legend.position='bottom',
      legend.title=element_blank(),
      strip.background= element_blank(),
      panel.border =    element_blank(),
      axis.title.x=     element_blank(),
      strip.text=       element_blank(),
      axis.title.y=     element_text(size=13),
      axis.text.y=      element_text(size=13), 
      axis.text.x=      element_text(size=13, angle = 45),
      legend.text =     element_text(size=12))


nhanes_data$AgeGroup <- as.character(nhanes_data$AgeGroup)
p1 <- nhanes_data |>
  select(AgeGroup, year, URXOPM) |> melt() |> 
  group_by(year, AgeGroup) |> 
  summarise(gm = geoMean(value, na.rm=T)) |>
  add_column(Cohort = c(rep("99-00", 3),
      rep(c("01-02", "07-08", "09-10", "11-12", "13-14"), each=4),
      rep("15-16", 5))) |>
  mutate(AgeGroup = replace(AgeGroup, AgeGroup == "Over 65", "> 65")) |>
  mutate(AgeGroup = replace(AgeGroup, AgeGroup == "0 - 5", "< 6")) |>
  mutate(across(AgeGroup, as.factor)) |>
  mutate(AgeGroup=fct_relevel(AgeGroup, levels(X$Age))) |> 
  mutate(across(Cohort, as.factor)) |>
  mutate(Cohort=fct_relevel(Cohort, levels(X$Cohort))) |> 
  ggplot(aes(x=Cohort, y = gm, group = AgeGroup)) + 
  geom_point(size=3) +
  geom_line() +
  facet_wrap(~AgeGroup, ncol=5, scale="free_x") +
  ylab("Urine 3PBA (ng/ml)") +
  theme_bw() + 
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=12),
    panel.border = element_blank(),
    axis.title.y=element_text(size=12),
    axis.title.x=element_blank(),
    axis.text.y=element_text(size=12),
    axis.text.x=element_text(size=12, angle = 45))

png("fig6_trend.png", width = 3000, height = 1500, res = 230)
plot_grid(p1, p2, nrow = 2, rel_heights = c(1/3, 2/3))
dev.off()

#
# Refer to https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/

dendro_data_k <- function(hc, k) {
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  hcdata
}

set_labels_params <- function(nbLabels,
                              direction = c("tb", "bt", "lr", "rl"),
                              fan       = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}

plot_ggdendro <- function(hcdata,
                          direction   = c("lr", "rl", "tb", "bt"),
                          fan         = FALSE,
                          scale.color = NULL,
                          branch.size = 1,
                          label.size  = 4,
                          nudge.label = 0.01,
                          expand.y    = 0.1) {
  # 
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)
  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 linewidth         =  branch.size) 
  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }
  # labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle
  p <- p +
    geom_text(data        =  label(hcdata),
              aes(x       =  x,
                  y       =  y,
                  label   =  label,
                  angle   =  angle,
                  colour  =  factor(clust)),
              #vjust       =  labelParams$vjust,
              #hjust       =  labelParams$hjust,
              vjust = 2.5,
              hjust = 1.0,
              nudge_y     =  ymax * nudge.label,
              size        =  label.size,
              show.legend =  FALSE)
  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }
  ylim <- -round(ymax * expand.y, 1)
  p <- p + expand_limits(y = ylim)
  p
}

hca_cyf <- X |> filter(Approach == "PBK" & Compound == "Cyfluthrin") |>
  select(Age, Compound, Cohort, Prediction) |>
  acast(Age~Cohort, median)|> scale() |>  
  dist(method = "euclidean") |> hclust(method = "ward.D2")
hca_dlm <- X |> filter(Approach == "PBK" & Compound == "Deltamethrin") |>
  select(Age, Compound, Cohort, Prediction) |>
  acast(Age~Cohort, median)|> scale() |>  
  dist(method = "euclidean") |> hclust(method = "ward.D2")
hca_cpm <- X |> filter(Approach == "PBK" & Compound == "Cypermethrin") |>
  select(Age, Compound, Cohort, Prediction) |>
  acast(Age~Cohort, median)|> scale() |>  
  dist(method = "euclidean") |> hclust(method = "ward.D2")
hca_prm <- X |> filter(Approach == "PBK" & Compound == "Permethrin") |>
  select(Age, Compound, Cohort, Prediction) |>
  acast(Age~Cohort, median)|> scale() |>  
  dist(method = "euclidean") |> hclust(method = "ward.D2")
hc_a <- X |> filter(Approach == "PBK") |>
  select(Age, Compound, Cohort, Prediction) |>
  unite("CC", Compound:Cohort, remove = T) |> 
  acast(Age~CC, median) |> scale() |>  
  dist(method = "euclidean") |> hclust(method = "ward.D2")

set_theme <- theme(plot.title = element_text(face='bold')) 

p1 <- dendro_data_k(hca_cyf, 2) |>  
  plot_ggdendro(direction = "tb", expand.y = 0.2) +
  theme_void() + ggtitle("Cyfluthrin") + set_theme
p2 <- dendro_data_k(hca_cpm, 2) |>  
  plot_ggdendro(direction = "tb", expand.y = 0.2) +
  theme_void() + ggtitle("Cypermethrin") + set_theme
p3 <- dendro_data_k(hca_dlm, 2) |>  
  plot_ggdendro(direction = "tb", expand.y = 0.2) +
  theme_void() + ggtitle("Deltamethrin") + set_theme
p4 <- dendro_data_k(hca_prm, 2) |>  
  plot_ggdendro(direction = "tb", expand.y = 0.2) +
  theme_void() + ggtitle("Permethrin") + set_theme
p5 <- dendro_data_k(hc_a, 2) |>  
  plot_ggdendro(direction = "tb", expand.y = 0.2) +
  theme_void() + ggtitle("Total") + set_theme

png("fig7_cluster.png", width = 2800, height = 2100, res = 300)
plot_grid(p5, plot_grid(p1, p2, p3, p4,  nrow = 1), nrow = 2)
dev.off()

