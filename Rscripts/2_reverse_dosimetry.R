# Load required libraries
library(readxl)     # read_excel
library(foreign)    # read.xport
library(survey)     # svydesign
library(magrittr)   # %<>%
library(dplyr)      # mutate
library(RMCSim)     # Monte Carlo simulations
library(foreach)    # Parallel loops
library(doParallel) # Parallel backend

# Load preprocessed data
load("data/Wtns.rda")
load("data/Agens.rda")
load("data/modparms.rda")

# Read NHANES metadata
codes_file <- "data/NHANEScodes_file.xlsx"
wtvars <- read_excel(codes_file, sheet = 2) |> as.data.frame()
demofiles <- setNames(wtvars$demofile, wtvars$sample)
datafiles <- setNames(wtvars$file, wtvars$sample)
bwtfiles <- setNames(wtvars$BWfile, wtvars$sample)
creatfiles <- setNames(wtvars$creatfile, wtvars$sample)
wtvar <- setNames(wtvars$wtvariable, wtvars$sample)

# Define cohorts and corresponding years
cohort <- c("99-00", "01-02", "07-08", "09-10", "11-12", "13-14", "15-16")
years <- c("1999", "2001", "2007", "2009", "2011", "2013", "2015")

# Copy executable for simulations
file.copy("MCSim/mcsim.gPYR_analytic_ss.model.exe", "mcsim.gPYR_analytic_ss.model.exe", overwrite = TRUE)

## Remove output folder or files ######
#system("rm -rf outputs")
#system("rm outputs/*.out")
#######################################

for (j in 1:7){
  
  select_cohort <- cohort[j]
  locat <- paste0("data/nhanes/", years[j])
  
  # File paths
  demo_locat <- paste0(
    locat, "/", unique(demofiles[which(names(demofiles) == cohort[j])])
  )
  datafiles <- c(
    "lab26pp.xpt", "l26pp_b.xpt", "uphopm_e.xpt", 
    "uphopm_f.xpt", "uphopm_g.xpt", "uphopm_h.xpt", "uphopm_i.xpt")
  data_locat <- paste0(locat, "/", datafiles[j])
  bwt_locat <- paste0(
    locat, "/", 
    tolower(unique(bwtfiles[which(names(bwtfiles) == cohort[j])]))
  )
  creat_locat <- paste0(
    locat, "/", unique(creatfiles[which(names(creatfiles) == cohort[j])]))
  
  # Read data
  demo <- read.xport(demo_locat)
  cdta <- read.xport(data_locat)
  bwt <- read.xport(bwt_locat)
  creat <- read.xport(creat_locat)
  
  # Define weight variable
  wtvariable  <- c(
    "WTSPP2YR", "WTSPP2YR", "WTSC2YR", "WTSC2YR",
    "WTSC2YR", "WTSC2YR", "WTSB2YR")
  chem2yrwt <- wtvariable[j]
  
  # Select chemical variables based on availability
  met <- c("URX4FP", "URXOPM", "URXTCC", "URXCB3", "URXCCC")
  if (all(met %in% names(cdta) == TRUE)) {
    chemvars <- c("URX4FP", "URXCB3", "URXCCC", "URXOPM", "URXTCC")
    LODnames <- c("URD4FPLC", "URDCB3LC", "URDCCCLC", "URDOPMLC", "URDTCCLC")
  } else if (all(met[1:4] %in% names(cdta) == TRUE)){
    chemvars <- c("URX4FP", "URXCB3", "URXOPM", "URXTCC")
    LODnames <- c("URD4FPLC", "URDCB3LC", "URDOPMLC", "URDTCCLC")
  } else {
    chemvars <- c("URX4FP", "URXOPM", "URXTCC") 
    LODnames <- c("URD4FPLC", "URDOPMLC", "URDTCCLC")
  }
  
  ## -----------------------------------------------------------------
  ## Demographics processing
  
  seq <- "SEQN"
  PSU <- "SDMVPSU" # the sampling unit
  STRA <- "SDMVSTRA" # the stratum for each observation
  demoageyr <- "RIDAGEYR"
  demogendr <- "RIAGENDR"
  demoeth <- "RIDRETH1"
  MECwt <- "WTMEC2YR"
  bodywtcomment <- "BMIWT"
    
  demo <- demo[,c(seq,PSU,STRA,demoageyr,demogendr,demoeth,MECwt)]
  ## Set up gender, age, and ethnicity as factors 
  # using the same levels as the NHANES reports
  demo[,demogendr] <- factor(demo[,demogendr], labels=c("Male","Female"))
  demo$AgeGroup <- cut(
    demo[,demoageyr], breaks=c(-1,5.5,11.5,19.5,65.5, 100.5),
    labels=c(
      "0 - 5 years", "6 - 11 years","12 - 19 years", 
      "20 - 65 years", "over 65 years"
    )
  )
  demo$RaceEthn <- factor(
    demo[,demoeth],
    labels=c(
      "Mexican American",
      "Other Hispanic","Non-Hispanic White",
      "Non-Hispanic Black","Other")
  )
  demo$ChildBearingAgeFemale <- factor(
    demo[,demogendr] == "Female" & (demo[,demoageyr] >= 16 & demo[,demoageyr] <= 49),
    labels=c("NotReproAgeFemale","ReproAgeFemale")
  )
  
  ## -----------------------------------------------------------------
  ## Chemical data processing  
  
  measurehead="URX"
  measuretail=NULL
  lodindhd="URD"
  lodtail="LC"
  
  meascore <- sub(
    paste("^",measurehead,"(.+)",measuretail,"$",sep=""),
    "\\1",
    chemvars
  )
  LODnames <- paste(lodindhd,meascore,lodtail,sep="")
  ## nms is a data frame with a record for each chemical.  It keeps
  ## track of the names of the measurement and LODind variables, the
  ## LOD, and space for the chemical name and CAS. LOD here is the
  ## maximum of the individual level LODs.
  nms <- data.frame(Measurement=chemvars,
                    LODind=LODnames,
                    LOD=numeric(length(chemvars)),
                    Chem = character(length(chemvars)),
                    CAS = character(length(chemvars)),
                    stringsAsFactors=FALSE)
  ## loop through the variables and find the lods.
  ## Need to confirm that there is only one value associated with an lod,
  ## then, the next to lowest value in the variable is the lod.
  for (i in seq_len(nrow(nms))) {
    nm <- nms$Measurement[i]
    nmlc <- nms$LODind[i]
    ## Do we already have an LOD indicator?
    if (!nmlc %in% colnames(cdta)) {
      ## We have to create the LOD indicator.  In this case, there is only 1 LOD
      ## and it will be sqrt(2) times the smallest value
      zlod <- sort(unique(cdta[,nm]))[1] * sqrt(2)
      cdta[,LODnames[i]] <- ifelse(cdta[,nm] >= zlod, 0, 1)
    }
  }
  
  #codes_file <- "data/NHANEScodes_file.xlsx"
  creatinine <- "URXUCR"
  
  ## -----------------------------------------------------------------
  ## Merge datasets
  alldata <- merge(demo, cdta[,c(seq, chem2yrwt, chemvars, LODnames)],
                   by.x=seq, by.y=seq, all.y=TRUE)
  
  alldata <- merge(alldata, creat[,c(seq, creatinine)],
                   by.x=seq, by.y=seq, all.x=TRUE) # 2125
  
  
  bodywt <- "BMXWT"
  bodymassindex <- "BMXBMI"
  
  
  bwt[!is.na(bwt[,bodywtcomment]),bodywt] <- NA
  ## Create obesity factor from bodymassindex
  bwt$Obesity <- cut(bwt[,bodymassindex], breaks=c(-0.5,30,500),
                     labels=c("BMI lteq 30", "BMI gt 30"))
  alldata <- merge(alldata, bwt[,c(seq, bodywt, bodymassindex, "Obesity")],
                   by.x=seq, by.y=seq, all.x=TRUE)
  
  ## -----------------------------------------------------------------
  ## Handle missing body weight
  if (any(is.na(alldata[,bodywt]))) {
    ##writeLines(paste(sum(is.na(alldata[,bodywt])),"missing values in",bodywt))
    ##browser()
    dta <- merge(demo, 
                 bwt[,c(seq, bodywt)], by.x=seq, by.y=seq, all.x=TRUE, all.y=TRUE
    )
    dsg <- na.omit(
      svydesign(ids=make.formula(PSU), strata=make.formula(STRA),
                weights=make.formula(MECwt), nest=TRUE, data=dta)
    ) 
    ## Just fit separate means by age and gender
    selectmales <- dsg$variables[,demogendr] == "Male"
    selectfemales <- dsg$variables[,demogendr] == "Female"
    males <- svyby(make.formula(bodywt), by=make.formula(demoageyr),
                   subset(dsg, selectmales), svymean)
    females <- svyby(make.formula(bodywt), by=make.formula(demoageyr),
                     subset(dsg, selectfemales), svymean)
    imp <- data.matrix(cbind(males[,bodywt],females[,bodywt]))
    ## Now, do the imputation
    ##browser()
    isnabw <- is.na(alldata[,bodywt])
    alldata[isnabw,bodywt] <- imp[data.matrix(alldata[isnabw,c(demoageyr, demogendr)])]
    ##writeLines(">>> Imputed bodyweights")
    ##print(alldata[isnabw, c(demoageyr, demogendr, bodywt)])
  }

  ## Fixup factors so there are no missing levels
  for (nm in names(alldata)) {
    if (is.factor(alldata[,nm])) {
      alldata[,nm] <- factor(alldata[,nm])
    }
  }
  
  ## -----------------------------------------------------------------
  # Calculate daily creatinine
  CreatFun <- function(newdata){
    X <- cbind(model.matrix(~ 0 + RIAGENDR + RIDRETH1, data=newdata),
               predict(Wtns, newdata$BMXWT),
               predict(Agens, newdata$RIDAGEYR))
    return(10^(X %*% modparms[-length(modparms)]))
  }
  
  alldata$DailyCreatinine <-
    CreatFun(data.frame(RIAGENDR = unname(alldata[,demogendr]),
                        RIDRETH1 = unname(alldata[,"RaceEthn"]),
                        BMXWT = unname(alldata[,bodywt]),
                        RIDAGEYR = unname(alldata[,demoageyr])))
  
  
  if (length(which(is.na(alldata$URXUCR))) > 0) 
    alldata <- alldata[-which(is.na(alldata$URXUCR)), ]
  if (length(which(is.na(alldata$DailyCreatinine))) > 0) 
    alldata <- alldata[-which(is.na(alldata$DailyCreatinine)), ]
  if (length(which(is.na(alldata$BMXWT))) > 0) 
    alldata <- alldata[-which(is.na(alldata$BMXWT)), ]
  
  ## -----------------------------------------------------------------
  # Handle non-detects
  LD_URX4FP <- min(alldata$URX4FP, na.rm = T)
  LD_URXTCC <- min(alldata$URXTCC, na.rm = T)
  LD_URXOPM <- min(alldata$URXOPM, na.rm = T)
  if ("URXCB3" %in% names(alldata)) 
    LD_URXCB3 <- min(alldata$URXCB3, na.rm = T)
  if ("URXCCC" %in% names(alldata)) 
    LD_URXCCC <- min(alldata$URXCCC, na.rm = T)
  
  ## Set NA to -1
  alldata$URX4FP[which(is.na(alldata$URX4FP))] <- -1
  alldata$URXTCC[which(is.na(alldata$URXTCC))] <- -1
  alldata$URXOPM[which(is.na(alldata$URXOPM))] <- -1
  if ("URXCB3" %in% names(alldata)) 
    alldata$URXCB3[which(is.na(alldata$URXCB3))] <- -1
  if ("URXCCC" %in% names(alldata)) 
    alldata$URXCCC[which(is.na(alldata$URXCCC))] <- -1
  
  #length(which(alldata$URX4FP == LD_URX4FP)) / length(alldata$URX4FP)
  #length(which(alldata$URXTCC == LD_URXTCC)) / length(alldata$URXTCC)
  #length(which(alldata$URXOPM == LD_URXOPM)) / length(alldata$URXOPM)
  #if ("URXCB3" %in% names(alldata)) 
  #  length(which(alldata$URXCB3 == LD_URXCB3)) / length(alldata$URXCB3)
  #if ("URXCCC" %in% names(alldata)) 
  #  length(which(alldata$URXCCC == LD_URXCCC)) / length(alldata$URXCCC)
  
  #
  alldata %<>% mutate(URX4FPND = ifelse(URX4FP == LD_URX4FP, LD_URX4FP*2^0.5/2, -1)) # set LOD = LOD/2 for t1
  #alldata %<>% mutate(URX4FPND = ifelse(URX4FP == LD_URX4FP,  0, -1)) # set LOD = 0
  alldata$URX4FP[which(alldata$URX4FP == LD_URX4FP)] <- -1
  alldata %<>% mutate(URXTCCND = ifelse(URXTCC == LD_URXTCC, LD_URXTCC*2^0.5/2, -1))
  #alldata %<>% mutate(URXTCCND = ifelse(URXTCC == LD_URXTCC, 0, -1))
  alldata$URXTCC[which(alldata$URXTCC == LD_URXTCC)] <- -1
  alldata %<>% mutate(URXOPMND = ifelse(URXOPM == LD_URXOPM, LD_URXOPM*2^0.5/2, -1))
  #alldata %<>% mutate(URXOPMND = ifelse(URXOPM == LD_URXOPM, 0, -1))
  alldata$URXOPM[which(alldata$URXOPM == LD_URXOPM)] <- -1
  if ("URXCB3" %in% names(alldata)){
    #alldata %<>% mutate(URXCB3ND = ifelse(URXCB3 == LD_URXCB3, 0.0, URXCB3))
    alldata %<>% mutate(URXCB3ND = ifelse(URXCB3 == LD_URXCB3, 0, -1))
    alldata$URXCB3[which(alldata$URXCB3 == LD_URXCB3)] <- -1
  } else {
    alldata$URXCB3ND <- 0.05 # set to prevent unrealistic prediction for t3
    alldata$URXCB3 <- -1
  }
  if ("URXCCC" %in% names(alldata)){
    #alldata %<>% mutate(URXCCCND = ifelse(URXCCC == LD_URXCCC, 0.0, URXCCC))
    alldata %<>% mutate(URXCCCND = ifelse(URXCCC == LD_URXCCC, 0, -1))
    alldata$URXCCC[which(alldata$URXCCC == LD_URXCCC)] <- -1
  } else {
    alldata$URXCCCND <- -1
    alldata$URXCCC <- -1
  }
 
  #alldata %<>% mutate(URX4FPND = ifelse(URX4FP == LD_URX4FP, (LD_URX4FP*2^0.5) / 2, -1)) # set LOD = LOD/2
  #alldata$URX4FP[which(alldata$URX4FP == LD_URX4FP)] <- -1
  #alldata %<>% mutate(URXTCCND = ifelse(URXTCC == LD_URXTCC, (LD_URXTCC*2^0.5) / 2, -1))
  #alldata$URXTCC[which(alldata$URXTCC == LD_URXTCC)] <- -1
  #alldata %<>% mutate(URXOPMND = ifelse(URXOPM == LD_URXOPM, (LD_URXOPM*2^0.5) / 2, -1))
  #alldata$URXOPM[which(alldata$URXOPM == LD_URXOPM)] <- -1
  #if ("URXCB3" %in% names(alldata)){
  #  #alldata %<>% mutate(URXCB3ND = ifelse(URXCB3 == LD_URXCB3, 0.0, URXCB3))
  #  alldata %<>% mutate(URXCB3ND = ifelse(URXCB3 == LD_URXCB3, (LD_URXCB3*2^0.5) / 2, -1))
  #  alldata$URXCB3[which(alldata$URXCB3 == LD_URXCB3)] <- -1
  #} else {
  #  alldata$URXCB3ND <- -1
  #  alldata$URXCB3 <- -1
  #}
  #if ("URXCCC" %in% names(alldata)){
  #  #alldata %<>% mutate(URXCCCND = ifelse(URXCCC == LD_URXCCC, 0.0, URXCCC))
  #  alldata %<>% mutate(URXCCCND = ifelse(URXCCC == LD_URXCCC, (LD_URXCCC*2^0.5) / 2, -1))
  #  alldata$URXCCC[which(alldata$URXCCC == LD_URXCCC)] <- -1
  #} else {
  #  alldata$URXCCCND <- -1
  #  alldata$URXCCC <- -1
  #}

  ## -----------------------------------------------------------------
  # Define groups for analysis
  total_groups <- length(unique(alldata$AgeGroup)) + 1

  gps <- c("Total", 
           "6 - 11 years", "12 - 19 years", "20 - 65 years",
           "over 65 years", "0 - 5 years",
           "BMI gt 30", "BMI lteq 30", 
           "Male", "Female", "ReproAgeFemale")
           #levels(alldata$AgeGroup)) 

  if (j %in% c(1:2)) {
    opt <- "MCSim/gPYR_analytic_5met.mcmc.in"
    opt_check <- "MCSim/gPYR_analytic_5met_check.mcmc.in"
  } else if (j %in% c(3:4)) {
    opt <- "MCSim/gPYR_analytic_4met.mcmc.in"
    opt_check <- "MCSim/gPYR_analytic_4met_check.mcmc.in"
  } else if (j %in% c(5:7)) {
    opt <- "MCSim/gPYR_analytic_3met.mcmc.in"
    opt_check <- "MCSim/gPYR_analytic_3met_check.mcmc.in"
  }
 
  #if (j %in% c(1:2)) {
  #  opt <- "MCSim/gPYR_analytic_5met.mcmc.in"
  #  opt_check <- "MCSim/gPYR_analytic_5met_check.mcmc.in"
  #} else if (j %in% c(3:4)) {
  #  opt <- "MCSim/gPYR_analytic_4met.mcmc.in"
  #  opt_check <- "MCSim/gPYR_analytic_4met_check.mcmc.in"
  #} else if (j %in% c(5:7)) {
  #  opt <- "MCSim/gPYR_analytic_3met.mcmc.in"
  #  opt_check <- "MCSim/gPYR_analytic_3met_check.mcmc.in"
  #}


  #for (k in c(1,7,8,9,10,11)){
  for (k in seq(total_groups)){
    if (k==1) {
      to_select <- c(1:dim(alldata)[1])
    } else if (k==2) {
      to_select <- which(alldata$AgeGroup == gps[2])
    } else if (k==3) {
      to_select <- which(alldata$AgeGroup == gps[3])
    } else if (k==4) {
      to_select <- which(alldata$AgeGroup == gps[4])
    } else if (k==5) {
      to_select <- which(alldata$AgeGroup == gps[5])
    } else if (k==6) {
      to_select <- which(alldata$AgeGroup == gps[6])
    } else if (k==7) {
      to_select <- which(alldata$Obesity == "BMI gt 30")
    } else if (k==8) {
      to_select <- which(alldata$Obesity == "BMI lteq 30")
    } else if (k==9) {
      to_select <- which(alldata$RIAGENDR == "Male")
    } else if (k==10) {
      to_select <- which(alldata$RIAGENDR == "Female")
    } else if (k==11) {
      to_select <- which(alldata$ChildBearingAgeFemale == "ReproAgeFemale")
    }
    
    
    if (j %in% c(1:2)) {
      ori_tx  <- readLines(opt)
      sink(opt, append = TRUE) 
      for(i in to_select){
        cat(paste0("Simulation {BW = ", alldata$BMXWT[i], "; ", 
                   "DailyCreatinine = ", alldata$DailyCreatinine[i], "; ",
                   "UrineCreatinine = ", alldata$URXUCR[i], "; ",
                   "Print(CUri3PBA_ND, 1); Data(CUri3PBA_ND, ", alldata$URXOPMND[i],  "); ",
                   "Print(CUri3PBA, 1); Data(CUri3PBA, ", alldata$URXOPM[i], "); ",
                   "Print(CUriFPBA_ND, 1); Data(CUriFPBA_ND, ", alldata$URX4FPND[i],  "); ",
                   "Print(CUriFPBA, 1); Data(CUriFPBA, ", alldata$URX4FP[i],  "); ",
                   "Print(CUriDBCA_ND, 1); Data(CUriDBCA_ND, ", alldata$URXCB3ND[i],  "); ",
                   "Print(CUriDBCA, 1); Data(CUriDBCA, ", alldata$URXCB3[i] ,");",
                   "Print(CUriDCCAc_ND, 1); Data(CUriDCCAc_ND, ", alldata$URXCCCND[i],  "); ",
                   "Print(CUriDCCAc, 1); Data(CUriDCCAc, ", alldata$URXCCC[i], "); ",
                   "Print(CUriDCCAt_ND, 1); Data(CUriDCCAt_ND, ", alldata$URXTCCND[i],  "); ",
                   "Print(CUriDCCAt, 1); Data(CUriDCCAt, ", alldata$URXTCC[i], "); }\n"
        )
        )
      }
      cat("}}\n") 
      sink()
      ori_chk_tx  <- readLines(opt_check)
      sink(opt_check, append = TRUE)
      #for(i in seq(nrow(alldata))){
      for(i in to_select){
        cat(paste0("Simulation {BW = ", alldata$BMXWT[i], "; ", 
                   "DailyCreatinine = ", alldata$DailyCreatinine[i], "; ",
                   "UrineCreatinine = ", alldata$URXUCR[i], "; ",
                   "Print(CUri3PBA_ND, 1); Data(CUri3PBA_ND, ", alldata$URXOPMND[i],  "); ",
                   "Print(CUri3PBA, 1); Data(CUri3PBA, ", alldata$URXOPM[i], "); ",
                   "Print(CUriFPBA_ND, 1); Data(CUriFPBA_ND, ", alldata$URX4FPND[i],  "); ",
                   "Print(CUriFPBA, 1); Data(CUriFPBA, ", alldata$URX4FP[i],  "); ",
                   "Print(CUriDBCA_ND, 1); Data(CUriDBCA_ND, ", alldata$URXCB3ND[i],  "); ",
                   "Print(CUriDBCA, 1); Data(CUriDBCA, ", alldata$URXCB3[i] ,");",
                   "Print(CUriDCCAc_ND, 1); Data(CUriDCCAc_ND, ", alldata$URXCCCND[i],  "); ",
                   "Print(CUriDCCAc, 1); Data(CUriDCCAc, ", alldata$URXCCC[i], "); ",
                   "Print(CUriDCCAt_ND, 1); Data(CUriDCCAt_ND, ", alldata$URXTCCND[i],  "); ",
                   "Print(CUriDCCAt, 1); Data(CUriDCCAt, ", alldata$URXTCC[i], "); }\n"
        )
        )
      }
      cat("}}\n")
      sink()
      #copy_file <- paste0("MCSim/gPYR_analytic_5met_", select_cohort, "_", gps[k], ".mcmc.in") 
      #copy_check_file <- paste0("MCSim/gPYR_analytic_5met_check_", select_cohort, "_", gps[k], ".mcmc.in")
      #invisible(suppressWarnings(file.remove(c(copy_file, copy_check_file))))
      #file.copy("MCSim/gPYR_analytic_5met.mcmc.in", copy_file)
      #file.copy("MCSim/gPYR_analytic_5met_check.mcmc.in", copy_check_file)
      #input <- paste0("gPYR_analytic_5met_", select_cohort, "_", gps[k], ".mcmc.in")
      copy_file <- paste0("MCSim/gPYR_analytic_5met_", select_cohort, "_", gps[k], ".mcmc.in") 
      copy_check_file <- paste0("MCSim/gPYR_analytic_5met_check_", select_cohort, "_", gps[k], ".mcmc.in")
      invisible(suppressWarnings(file.remove(c(copy_file, copy_check_file))))
      file.copy("MCSim/gPYR_analytic_5met.mcmc.in", copy_file)
      file.copy("MCSim/gPYR_analytic_5met_check.mcmc.in", copy_check_file)
      input <- paste0("gPYR_analytic_5met_", select_cohort, "_", gps[k], ".mcmc.in")
    } else if (j %in% c(3:4)) {
      ori_tx  <- readLines(opt)
      sink(opt, append = TRUE) 
      for(i in to_select){
        cat(paste0("Simulation {BW = ", alldata$BMXWT[i], "; ", 
                   "DailyCreatinine = ", alldata$DailyCreatinine[i], "; ",
                   "UrineCreatinine = ", alldata$URXUCR[i], "; ",
                   "Print(CUri3PBA_ND, 1); Data(CUri3PBA_ND, ", alldata$URXOPMND[i],  "); ",
                   "Print(CUri3PBA, 1); Data(CUri3PBA, ", alldata$URXOPM[i], "); ",
                   "Print(CUriFPBA_ND, 1); Data(CUriFPBA_ND, ", alldata$URX4FPND[i],  "); ",
                   "Print(CUriFPBA, 1); Data(CUriFPBA, ", alldata$URX4FP[i],  "); ",
                   "Print(CUriDBCA_ND, 1); Data(CUriDBCA_ND, ", alldata$URXCB3ND[i],  "); ",
                   "Print(CUriDBCA, 1); Data(CUriDBCA, ", alldata$URXCB3[i] ,");",
                   #"Print(CUriDCCAc_ND, 1); Data(CUriDCCAc_ND, ", alldata$URXCCCND[i],  "); ",
                   #"Print(CUriDCCAc, 1); Data(CUriDCCAc, ", alldata$URXCCC[i], "); ",
                   "Print(CUriDCCAt_ND, 1); Data(CUriDCCAt_ND, ", alldata$URXTCCND[i],  "); ",
                   "Print(CUriDCCAt, 1); Data(CUriDCCAt, ", alldata$URXTCC[i], "); }\n"
        )
        )
      }
      cat("}}\n") 
      sink()
      ori_chk_tx  <- readLines(opt_check)
      sink(opt_check, append = TRUE)
      #for(i in seq(nrow(alldata))){
      for(i in to_select){
        cat(paste0("Simulation {BW = ", alldata$BMXWT[i], "; ", 
                   "DailyCreatinine = ", alldata$DailyCreatinine[i], "; ",
                   "UrineCreatinine = ", alldata$URXUCR[i], "; ",
                   "Print(CUri3PBA_ND, 1); Data(CUri3PBA_ND, ", alldata$URXOPMND[i],  "); ",
                   "Print(CUri3PBA, 1); Data(CUri3PBA, ", alldata$URXOPM[i], "); ",
                   "Print(CUriFPBA_ND, 1); Data(CUriFPBA_ND, ", alldata$URX4FPND[i],  "); ",
                   "Print(CUriFPBA, 1); Data(CUriFPBA, ", alldata$URX4FP[i],  "); ",
                   "Print(CUriDBCA_ND, 1); Data(CUriDBCA_ND, ", alldata$URXCB3ND[i],  "); ",
                   "Print(CUriDBCA, 1); Data(CUriDBCA, ", alldata$URXCB3[i] ,");",
                   #"Print(CUriDCCAc_ND, 1); Data(CUriDCCAc_ND, ", alldata$URXCCCND[i],  "); ",
                   #"Print(CUriDCCAc, 1); Data(CUriDCCAc, ", alldata$URXCCC[i], "); ",
                   "Print(CUriDCCAt_ND, 1); Data(CUriDCCAt_ND, ", alldata$URXTCCND[i],  "); ",
                   "Print(CUriDCCAt, 1); Data(CUriDCCAt, ", alldata$URXTCC[i], "); }\n"
        )
        )
      }
      cat("}}\n")
      sink()
      #copy_file <- paste0("MCSim/gPYR_analytic_4met_t1_", select_cohort, "_", gps[k], ".mcmc.in") 
      #copy_check_file <- paste0("MCSim/gPYR_analytic_4met_t1_check_", select_cohort, "_", gps[k], ".mcmc.in")
      #invisible(suppressWarnings(file.remove(c(copy_file, copy_check_file))))
      #file.copy("MCSim/gPYR_analytic_4met_t1.mcmc.in", copy_file)
      #file.copy("MCSim/gPYR_analytic_4met_t1_check.mcmc.in", copy_check_file)
      #input <- paste0("gPYR_analytic_4met_t1_", select_cohort, "_", gps[k], ".mcmc.in")
      copy_file <- paste0("MCSim/gPYR_analytic_4met_", select_cohort, "_", gps[k], ".mcmc.in") 
      copy_check_file <- paste0("MCSim/gPYR_analytic_4met_check_", select_cohort, "_", gps[k], ".mcmc.in")
      invisible(suppressWarnings(file.remove(c(copy_file, copy_check_file))))
      file.copy("MCSim/gPYR_analytic_4met.mcmc.in", copy_file)
      file.copy("MCSim/gPYR_analytic_4met_check.mcmc.in", copy_check_file)
      input <- paste0("gPYR_analytic_4met_", select_cohort, "_", gps[k], ".mcmc.in")
    } else if (j %in% c(5:7)) {
      ori_tx  <- readLines(opt)
      sink(opt, append = TRUE)
      for(i in to_select){
        cat(paste0("Simulation {BW = ", alldata$BMXWT[i], "; ", 
                   "DailyCreatinine = ", alldata$DailyCreatinine[i], "; ",
                   "UrineCreatinine = ", alldata$URXUCR[i], "; ",
                   "Print(CUri3PBA_ND, 1); Data(CUri3PBA_ND, ", alldata$URXOPMND[i],  "); ",
                   "Print(CUri3PBA, 1); Data(CUri3PBA, ", alldata$URXOPM[i], "); ",
                   "Print(CUriFPBA_ND, 1); Data(CUriFPBA_ND, ", alldata$URX4FPND[i],  "); ",
                   "Print(CUriFPBA, 1); Data(CUriFPBA, ", alldata$URX4FP[i],  "); ",
                   "Print(CUriDBCA_ND, 1); Data(CUriDBCA_ND, ", alldata$URXCB3ND[i],  "); ",
                   #"Print(CUriDBCA, 1); Data(CUriDBCA, ", alldata$URXCB3[i] ,");",
                   #"Print(CUriDCCAc_ND, 1); Data(CUriDCCAc_ND, ", alldata$URXCCCND[i],  "); ",
                   #"Print(CUriDCCAc, 1); Data(CUriDCCAc, ", alldata$URXCCC[i], "); ",
                   "Print(CUriDCCAt_ND, 1); Data(CUriDCCAt_ND, ", alldata$URXTCCND[i],  "); ",
                   "Print(CUriDCCAt, 1); Data(CUriDCCAt, ", alldata$URXTCC[i], "); }\n"
        )
        )
      }
      cat("}}\n") 
      sink()
      ori_chk_tx  <- readLines(opt_check)
      sink(opt_check, append = TRUE)
      #for(i in seq(nrow(alldata))){
      for(i in to_select){
        cat(paste0("Simulation {BW = ", alldata$BMXWT[i], "; ", 
                   "DailyCreatinine = ", alldata$DailyCreatinine[i], "; ",
                   "UrineCreatinine = ", alldata$URXUCR[i], "; ",
                   "Print(CUri3PBA_ND, 1); Data(CUri3PBA_ND, ", alldata$URXOPMND[i],  "); ",
                   "Print(CUri3PBA, 1); Data(CUri3PBA, ", alldata$URXOPM[i], "); ",
                   "Print(CUriFPBA_ND, 1); Data(CUriFPBA_ND, ", alldata$URX4FPND[i],  "); ",
                   "Print(CUriFPBA, 1); Data(CUriFPBA, ", alldata$URX4FP[i],  "); ",
                   "Print(CUriDBCA_ND, 1); Data(CUriDBCA_ND, ", alldata$URXCB3ND[i],  "); ",
                   #"Print(CUriDBCA, 1); Data(CUriDBCA, ", alldata$URXCB3[i] ,");",
                   #"Print(CUriDCCAc_ND, 1); Data(CUriDCCAc_ND, ", alldata$URXCCCND[i],  "); ",
                   #"Print(CUriDCCAc, 1); Data(CUriDCCAc, ", alldata$URXCCC[i], "); ",
                   "Print(CUriDCCAt_ND, 1); Data(CUriDCCAt_ND, ", alldata$URXTCCND[i],  "); ",
                   "Print(CUriDCCAt, 1); Data(CUriDCCAt, ", alldata$URXTCC[i], "); }\n"
        )
        )
      }
      cat("}}\n")
      sink()
      #copy_file <- paste0("MCSim/gPYR_analytic_3met_t1_", select_cohort, "_", gps[k], ".mcmc.in") 
      #copy_check_file <- paste0("MCSim/gPYR_analytic_3met_t1_check_", select_cohort, "_", gps[k], ".mcmc.in")
      #invisible(suppressWarnings(file.remove(c(copy_file, copy_check_file))))
      #file.copy("MCSim/gPYR_analytic_3met_t1.mcmc.in", copy_file)
      #file.copy("MCSim/gPYR_analytic_3met_t1_check.mcmc.in", copy_check_file)
      #input <- paste0("gPYR_analytic_3met_t1_", select_cohort, "_", gps[k], ".mcmc.in")
      copy_file <- paste0("MCSim/gPYR_analytic_3met_", select_cohort, "_", gps[k], ".mcmc.in") 
      copy_check_file <- paste0("MCSim/gPYR_analytic_3met_check_", select_cohort, "_", gps[k], ".mcmc.in")
      invisible(suppressWarnings(file.remove(c(copy_file, copy_check_file))))
      file.copy("MCSim/gPYR_analytic_3met.mcmc.in", copy_file)
      file.copy("MCSim/gPYR_analytic_3met_check.mcmc.in", copy_check_file)
      input <- paste0("gPYR_analytic_3met_", select_cohort, "_", gps[k], ".mcmc.in")
    }
   
    writeLines(ori_tx, con = opt)
    writeLines(ori_chk_tx, con = opt_check)
    
    
    
    #cat("The total number of participants in", select_cohort, ": ", length(to_select), "\n")
    # source("3_reverse_dosimetry")
    
    
    model <- "gPYR_analytic_ss.model"
    
    
    #system("rm -rf outputs")
    if (!dir.exists("outputs")) dir.create("outputs")
    
    current_files <- list.files()
    detectCores()
    cores <- 3    
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    strt <- Sys.time()
    out <- foreach(i = 1:cores) %dopar% {
      set.seed(i + 1)
      RMCSim::mcsim(model = model, input = input, dir = "MCSim", parallel = T, check=F)
    }
    print(Sys.time() - strt)
    
    new_files <- setdiff(list.files(), current_files)
    to_remove <- new_files[grep(".kernel|.in", new_files)]
    file.remove(to_remove)
    out_files <- setdiff(list.files(), current_files)
    
    for (i in 1:cores) file.copy(out_files[i], paste0("outputs/", out_files[i]))
    #file.remove(out_files) 
    system("rm *.out")

  }
}

file.remove("mcsim.gPYR_analytic_ss.model.exe")




