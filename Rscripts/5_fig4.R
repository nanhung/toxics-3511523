library(dplyr) # bind_rows
all_cohorts_data <- c("99-00", "01-02", 
                      #"03-04", "05-06", 
                      "07-08", "09-10", 
                      "11-12", "13-14", "15-16")


pyr_P <-c("68359-37-5", # CYF
          "52645-53-1", # PRM
          "52918-63-5", # DLM
          "52315-07-8") # CPM


#
#cohort <- "99-00"
for (cohort in all_cohorts_data){
  #
  OnlyPparms3 <- list.files(paste0("./saves/", cohort,"/"), pattern='lPsamps-gm')
  OnlyPparms3
  PBPK_files <- list.files(paste0("./saves/", cohort,"/"), pattern='PBPK')
  PBPK_files
  #
  op_gp_no <- c(
    grep("20 - 65", x = OnlyPparms3), 
    grep("12 - 19", x = OnlyPparms3), 
    grep("6 - 11", x = OnlyPparms3), 
    grep("66 years", x = OnlyPparms3), 
    grep("0 - 5", x = OnlyPparms3))
  pk_gp_no <- c(
    grep("20 - 65", x = PBPK_files), 
    grep("12 - 19", x = PBPK_files), 
    grep("6 - 11", x = PBPK_files),
    grep("over 65", x = PBPK_files), 
    grep("0 - 5", x = PBPK_files))
  op_gps <- OnlyPparms3[op_gp_no]
  op_gps
  pk_gps <- PBPK_files[pk_gp_no]
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
    OnlyPparms3_to_be_load <- paste0("saves/", cohort, "/", op_gps[i])
    OnlyPparms3_to_be_load
    load(OnlyPparms3_to_be_load)
    lPsampsgm |> dim()
    lPsampsgm <- lPsampsgm[randsamp, ]
    #
    PBPK_file_to_be_load <- paste0("saves/", cohort, "/", pk_gps[i])
    PBPK_file_to_be_load
    load(PBPK_file_to_be_load)
    lPsampsPBPK |> dim()
    lPsampsPBPK <- lPsampsPBPK[randsamp, ]
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
                           "Cyfluthrin", "High-throughput TK", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app1_prm <- data.frame(c(exp(lPsampsgm[,which(colnames(lPsampsgm) == pyr_P[2])])),
                           "Permethrin", "High-throughput TK", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app1_dlm <- data.frame(c(exp(lPsampsgm[,which(colnames(lPsampsgm) == pyr_P[3])])),
                           "Deltamethrin", "High-throughput TK", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app1_cpm <- data.frame(c(exp(lPsampsgm[,which(colnames(lPsampsgm) == pyr_P[4])])),
                           "Cypermethrin", "High-throughput TK", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app2_cyf <- data.frame(c(exp(lPsampsPBPK[,"lnCYF_IngDose(1)"])),
                           "Cyfluthrin", "PBK", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app2_prm <- data.frame(c(exp(lPsampsPBPK[,"lnPRM_IngDose(1)"])),
                           "Permethrin", "PBK", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    app2_cpm <- data.frame(c(exp(lPsampsPBPK[,"lnPRM_IngDose(1)"])*0.1),
                           "Cypermethrin", "PBK", cohort, age) |>
      `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    #if(cohort %in% all_cohorts_data[c(1:4)]){
      app2_dlm <- data.frame(c(exp(lPsampsPBPK[,"lnDLM_IngDose(1)"])),
                             "Deltamethrin", "PBK", cohort, age) |>
        `colnames<-`(c("Prediction", "Compound", "Approach", "Cohort", "Age"))
    #}
    #
    pbpk_pop_cyf <- data.frame(
      c(lPsampsPBPK[,"lnCYF_IngDose(1)"]),
      c(lPsampsPBPK[,"V_lnCYF_IngDose(1)"]),
      "Cyfluthrin", cohort, age) |>
      `colnames<-`(c("mean", "variance", "Compound", "Cohort", "Age"))
    pbpk_pop_prm <- data.frame(
      c(lPsampsPBPK[,"lnPRM_IngDose(1)"]),
      c(lPsampsPBPK[,"V_lnPRM_IngDose(1)"]),
      "Permethrin", cohort, age) |>
      `colnames<-`(c("mean", "variance", "Compound", "Cohort", "Age"))
    #if(cohort %in% all_cohorts_data[c(1:4)]){
      pbpk_pop_dlm <- data.frame(
        c(lPsampsPBPK[,"lnDLM_IngDose(1)"]),
        c(lPsampsPBPK[,"V_lnDLM_IngDose(1)"]),
        "Deltamethrin", cohort, age) |>
      `colnames<-`(c("mean", "variance", "Compound", "Cohort", "Age"))
    #}
    #
    if ( cohort=="99-00" & i==1 ){
      X <- bind_rows(app1_cyf, app2_cyf, 
                     app1_dlm, app2_dlm,
                     app1_prm, app2_prm,
                     app1_cpm, app2_cpm)
      X_pbpk <- bind_rows(pbpk_pop_cyf, pbpk_pop_prm, pbpk_pop_dlm)
    } else{
      X <- bind_rows(X,
                     app1_cyf, app2_cyf, 
                     app1_dlm, app2_dlm,
                     app1_prm, app2_prm,
                     app1_cpm, app2_cpm)
      X_pbpk <- bind_rows(X_pbpk, pbpk_pop_cyf, pbpk_pop_prm, pbpk_pop_dlm)
    }
  }
}
