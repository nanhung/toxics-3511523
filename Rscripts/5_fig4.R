library(dplyr) # bind_rows
library(ggpubr) # ggscatter
library(scales)

all_cohorts_data <- c("99-00", "01-02", 
                      #"03-04", "05-06", 
                      "07-08", "09-10", 
                      "11-12", "13-14", "15-16")


pyr_P <-c("68359-37-5", # CYF
          "52645-53-1", # PRM
          "52918-63-5", # DLM
          "52315-07-8") # CPM


#
#cohort <- "11-12"
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

png("fig4_comparison.png", width = 2000, height = 1500, res = 300)
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

