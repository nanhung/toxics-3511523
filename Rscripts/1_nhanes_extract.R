# The function was modified by bayesmarker package
# The NHANES recently updated the data site (02/2025) 

library(readxl) # read_excel
library(foreign) # read.xport

get_NHANES_data <- function(codes_file = NULL, cohort = NULL, save_directory = NULL) {
  if (is.null(codes_file)){
    print("Error: please provide a file name")
    stop()
  }
  codes <- as.data.frame(read_excel(codes_file, sheet = 1))
  weights <- as.data.frame(read_excel(codes_file, sheet = 2))
  # Generate list of needed laboratory files 
  # (has metabolite concentration measurements)
  samps <- unique(codes$recent_sample)
  # Limit by cohort
  if (!is.null(cohort)){
    samps <- samps[samps %in% cohort]
  }
  labFiles <- list()
  labFiles <- lapply(samps, function(x) 
    unique(codes$NHANESfile[codes$recent_sample == x]))
  names(labFiles) <- samps
  # Reference table for phases and file conventions
  ref <- data.frame("Full" = c("1999", "2001", "2003", 
                               "2005", "2007", "2009",
                               "2011", "2013", "2015", 
                               "2017", "2019"),
                    "Short" = c("99-00", "01-02", "03-04", "05-06", "07-08", 
                                "09-10", "11-12", "13-14", "15-16", "17-18", 
                                "19-20"),
                    "Ext" = c("", "B", "C", "D", "E", "F", "G", "H",
                              "I", "J", "K"),
                    stringsAsFactors = FALSE)
  ### Add demographic, bodyweight, and creatinine files
  result <- list()
  for (i in 1:length(samps)) {
    ind <- match(samps[i], weights$sample)
    result[[i]] <- c(paste(labFiles[[i]], ".XPT", sep = ""),
                     paste(ifelse(samps[i] != "99-00", "DEMO_", "DEMO"), 
                           ref$Ext[ref$Short == samps[i]], ".XPT", sep = ""),
                     paste(ifelse(samps[i] != "99-00", "BMX_", "BMX"), 
                           ref$Ext[ref$Short == samps[i]], ".XPT", sep = ""),
                     toupper(weights$creatfile[ind]), 
                     toupper(weights$urineflow[ind]))
  }
  names(result) <- samps 
  ### ??? debug to solve the issue for 99-00, 01-02, 03-04 data
  for (i in 1:length(result)) 
    result[[i]] <- result[[i]][complete.cases(result[[i]])]
  ### ???
  ### Download the needed files from the NHANES website and save in the correct directory
  print("Starting Downlodas")
  oldw <- getOption("warn")
  options(warn = -1)
  for (i in 1:length(samps)) {
    indP <- match(samps[i], ref$Short)
    # Create nhanes directory
    if (is.null(save_directory)){
      if (!dir.exists("nhanes")) {
        dir.create("./nhanes")
      }
    } else {
      if (!dir.exists(save_directory)){
        dir.create(save_directory)
      }
      if (!dir.exists(file.path(save_directory, "nhanes/"))) {
        dir.create(file.path(save_directory, "nhanes/"))
      }
    }
    for (j in 1:length(result[[i]])) {
        if (result[[i]][j] == ""){
          next
        }
        if (is.null(save_directory)){
          if (!dir.exists(paste("nhanes/", ref$Full[indP], sep = ""))) {
            dir.create(paste("./nhanes", ref$Full[indP], sep = ""))
          }
          print(
            paste("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", 
              ref$Full[indP], "/DataFiles/", result[[i]][j], sep = "")
            )
          tryCatch(download.file(
            paste("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", 
              ref$Full[indP], "/DataFiles/", result[[i]][j], sep = ""), 
            destfile = paste("./nhanes/", ref$Full[indP], "/", 
              tolower(result[[i]][j]), sep = ""), mode = "wb"), 
          error = function(e) print(paste(result[[i]][j], 'was not found',
              sep = " ")))
        } else {
          if (!dir.exists(file.path(save_directory, 
                paste("nhanes/", ref$Full[indP], sep = "")))) {
           dir.create(file.path(save_directory, 
               paste("nhanes/", ref$Full[indP], sep = "")))
          }
          print(
            paste("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", 
              ref$Full[indP], "/DataFiles/", result[[i]][j], sep = "")
            )
          tryCatch(download.file(
            paste("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", 
              ref$Full[indP], "/DataFiles/", result[[i]][j], sep = ""),
            destfile = file.path(save_directory,
              paste("nhanes/", ref$Full[indP], "/", tolower(result[[i]][j]), sep = "")),
            mode = "wb"), error = function(e) print(paste(result[[i]][j], 'was not found', sep = " ")))
        }
      }
  }
  options(warn = oldw)
}


get_data <- function(cohort = "99-00"){
  get_NHANES_data(codes_file = "data/NHANEScodes_file.xlsx",
                  cohort = cohort,
                  save_directory = "data")
}

years <- c("1999", "2001", "2007", "2009", "2011", "2013", "2015")
cohort <- c("99-00", "01-02", "07-08", "09-10", "11-12", "13-14", "15-16")

for(i in seq(7)) {
  check_files <- list.files(paste0("data/nhanes/", years[i], "/"))
  if (length(check_files) == 0) get_data(cohort = cohort[i])
}

#
codes_file <- "data/NHANEScodes_file.xlsx"
NHANEScodes <- codes_file
wtvars <- as.data.frame(read_excel(NHANEScodes, sheet = 2))
demofiles <- wtvars$demofile
datafiles <- wtvars$file
bwtfiles <- wtvars$BWfile
creatfiles <- wtvars$creatfile
wtvar <- wtvars$wtvariable
names(demofiles) <- names(datafiles) <-
  names(bwtfiles) <- names(creatfiles) <- 
  names(wtvar) <- wtvars$sample

for (j in seq(7)) {
  select_cohort <- cohort[j]
  locat <- paste0("data/nhanes/", years[j])
  demo_locat <- paste0(
    locat, "/", unique(demofiles[which(names(demofiles) == cohort[j])])
  )
  datafiles <- c("lab26pp.xpt", "l26pp_b.xpt", "uphopm_e.xpt", 
    "uphopm_f.xpt", "uphopm_g.xpt", "uphopm_h.xpt", "uphopm_i.xpt")
  data_locat <- paste0(locat, "/", datafiles[j])
  bwt_locat <- paste0(
    locat, "/", tolower(unique(bwtfiles[which(names(bwtfiles) == cohort[j])]))
  )
  creat_locat <- paste0(
    locat, "/",unique(creatfiles[which(names(creatfiles) == cohort[j])])
  )
  demo <- read.xport(demo_locat)
  cdta <- read.xport(data_locat)
  bwt <- read.xport(bwt_locat)
  creat <- read.xport(creat_locat)
  wtvariable  <- c("WTSPP2YR", "WTSPP2YR", "WTSC2YR", 
    "WTSC2YR", "WTSC2YR", "WTSC2YR", "WTSB2YR")
  chem2yrwt <- wtvariable[j]
  met <- c("URX4FP", "URXOPM", "URXTCC", "URXCB3", "URXCCC")
  if (all(met %in% names(cdta) == TRUE)) {
    chemvars <- c("URX4FP", "URXOPM", "URXTCC", "URXCB3", "URXCCC")
    LODnames <- c("URD4FPLC", "URDCB3LC", "URDCCCLC", "URDOPMLC", "URDTCCLC")
  } else if (all(met[1:4] %in% names(cdta) == TRUE)){
    chemvars <- c("URX4FP", "URXOPM", "URXTCC", "URXCB3")
    LODnames <- c("URD4FPLC", "URDCB3LC", "URDOPMLC", "URDTCCLC")
  } else {
    chemvars <- c("URX4FP", "URXOPM", "URXTCC") 
    LODnames <- c("URD4FPLC", "URDOPMLC", "URDTCCLC")
  }
  seq <- "SEQN"
  PSU <- "SDMVPSU" # the sampling unit
  STRA <- "SDMVSTRA" # the stratum for each observation
  demoageyr <- "RIDAGEYR"
  demogendr <- "RIAGENDR"
  demoeth <- "RIDRETH1"
  MECwt <- "WTMEC2YR"
  bodywtcomment <- "BMIWT"
  demo <- demo[,c(seq,PSU,STRA,demoageyr,demogendr,demoeth,MECwt)]
  # Set up gender, age, and ethnicity as factors 
  # using the same levels as the NHANES reports
  demo[,demogendr] <- factor(demo[,demogendr], labels=c("Male","Female"))
  demo$AgeGroup <- cut(
    demo[,demoageyr], 
    breaks=c(-1,5.5,11.5,19.5,65.5, 100.5),
    labels=c("0 - 5", "6 - 11", "12 - 19", "20 - 65", "Over 65")
  )
  demo$RaceEthn <- factor(
    demo[,demoeth],
    labels=c(
      "Mexican American",
      "Other Hispanic",
      "Non-Hispanic White",
      "Non-Hispanic Black","Other")
  )
  demo$ChildBearingAgeFemale <- factor(
    demo[,demogendr] == "Female" & (demo[,demoageyr] >= 16 & demo[,demoageyr] <= 49),
    labels=c("NotReproAgeFemale","ReproAgeFemale")
  )
  measurehead="URX"
  measuretail=NULL
  lodindhd="URD"
  lodtail="LC"
  meascore <- sub(
    paste("^",measurehead,"(.+)",measuretail,"$",sep=""),"\\1",chemvars
  )
  LODnames <- paste(lodindhd,meascore,lodtail,sep="")
  ## nms is a data frame with a record for each chemical.  It keeps
  ## track of the names of the measurement and LODind variables, the
  ## LOD, and space for the chemical name and CAS. LOD here is the
  ## maximum of the individual level LODs.
  nms <- data.frame(
    Measurement=chemvars, 
    LODind=LODnames, 
    LOD=numeric(length(chemvars)),
    Chem = character(length(chemvars)), 
    CAS = character(length(chemvars)),
    stringsAsFactors=FALSE
  )
  for (i in seq_len(nrow(nms))) {
    nm <- nms$Measurement[i]
    nmlc <- nms$LODind[i]
    ## Do we already have an LOD indicator?
    if (!nmlc %in% colnames(cdta)) {
      ## We have to create the LOD indicator.  
      ## In this case, there is only 1 LOD
      ## and it will be sqrt(2) times the smallest value
      zlod <- sort(unique(cdta[,nm]))[1] * sqrt(2)
      cdta[,LODnames[i]] <- ifelse(cdta[,nm] >= zlod, 0, 1)
    }
  }
  codes_file <- "data/NHANEScodes_file.xlsx"
  creatinine <- "URXUCR"
  alldata <- merge(demo, cdta[,c(seq, chem2yrwt, chemvars, LODnames)],
    by.x=seq, by.y=seq, all.y=TRUE)
  alldata <- merge(alldata, creat[,c(seq, creatinine)],
    by.x=seq, by.y=seq, all.x=TRUE) # 2125
  bodywt <- "BMXWT"
  bodymassindex <- "BMXBMI"
  bwt[!is.na(bwt[,bodywtcomment]),bodywt] <- NA
  ## Create obesity factor from bodymassindex
  bwt$Obesity <- cut(bwt[,bodymassindex], breaks=c(-0.5,30,500),
    labels=c("BMI <= 30", "BMI > 30"))
  alldata <- merge(alldata, bwt[,c(seq, bodywt, bodymassindex, "Obesity")],
    by.x=seq, by.y=seq, all.x=TRUE)
  if (any(is.na(alldata[,bodywt]))) {
    ##writeLines(paste(sum(is.na(alldata[,bodywt])),"missing values in",bodywt))
    ##browser()
    dta <- merge(
      demo, bwt[,c(seq, bodywt)], by.x=seq, by.y=seq, all.x=TRUE, all.y=TRUE)
    dsg <- na.omit(svydesign(ids=make.formula(PSU), strata=make.formula(STRA),
        weights=make.formula(MECwt), nest=TRUE, data=dta))
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
    alldata[isnabw,bodywt] <- imp[
      data.matrix(alldata[isnabw,c(demoageyr, demogendr)])
      ]
    ##writeLines(">>> Imputed bodyweights")
    ##print(alldata[isnabw, c(demoageyr, demogendr, bodywt)])
  }
  ## Fixup factors so there are no missing levels
  for (nm in names(alldata)) {
    if (is.factor(alldata[,nm])) {
      alldata[,nm] <- factor(alldata[,nm])}
  }
  #
  CreatFun <- function(newdata){
    X <- cbind(model.matrix(~ 0 + RIAGENDR + RIDRETH1, data=newdata),
      predict(Wtns, newdata$BMXWT),
      predict(Agens, newdata$RIDAGEYR))
    return(10^(X %*% modparms[-length(modparms)]))
  }
  #
  load("data/Wtns.rda")
  load("data/Agens.rda")
  load("data/modparms.rda")
  alldata$DailyCreatinine <-
    CreatFun(data.frame(RIAGENDR = unname(alldata[,demogendr]),
        RIDRETH1 = unname(alldata[,"RaceEthn"]),
        BMXWT = unname(alldata[,bodywt]),
        RIDAGEYR = unname(alldata[,demoageyr])))
  #
  LD_URXOPM <- min(alldata$URXOPM, na.rm = T)
  alldata$LD_URXOPM <- min(alldata$URXOPM, na.rm = T)
  LD_URXTCC <- min(alldata$URXTCC, na.rm = T)
  alldata$LD_URXTCC <- min(alldata$URXTCC, na.rm = T)
  LD_URX4FP <- min(alldata$URX4FP, na.rm = T)
  alldata$LD_URX4FP <- min(alldata$URX4FP, na.rm = T)
  if (all(met %in% names(cdta) == TRUE)){
    LD_URXCB3 <- min(alldata$URXCB3, na.rm = T)
    alldata$LD_URXCB3 <- min(alldata$URXCB3, na.rm = T)
    LD_URXCCC <- min(alldata$URXCCC, na.rm = T)
    alldata$LD_URXCCC <- min(alldata$URXCCC, na.rm = T)
  } else if (all(met[1:4] %in% names(cdta) == TRUE)) {
    LD_URXCB3 <- min(alldata$URXCB3, na.rm = T) 
    alldata$LD_URXCB3 <- min(alldata$URXCB3, na.rm = T)
  }
  if (j %in% c(3:7)){
    alldata$URXCCC <- NA
    alldata$LD_URXCCC <- NA
  }
  if (j %in% c(5:7)){
    alldata$URXCB3 <- NA
    alldata$LD_URXCB3 <- NA
  }
  # tidy
  if (length(which(is.na(alldata$URXUCR))) > 0) {
    alldata <- alldata[-which(is.na(alldata$URXUCR)), ]
  }
  if (length(which(is.na(alldata$DailyCreatinine))) > 0) {
    alldata <- alldata[-which(is.na(alldata$DailyCreatinine)), ]
  }
  if (length(which(is.na(alldata$BMXWT))) > 0) {
    alldata <- alldata[-which(is.na(alldata$BMXWT)), ]
  }
  alldata$year <- years[j]
  to_select <- c("AgeGroup", "BMXWT", "DailyCreatinine", "URXUCR", 
    "URX4FP", "URXOPM", "URXTCC", "URXCB3", "URXCCC",
    "LD_URX4FP", "LD_URXOPM", "LD_URXTCC", "LD_URXCB3",  "LD_URXCCC", "year")
  alldata <- alldata[, to_select] 
  if (j == 1) all_data <- alldata else all_data <- rbind(all_data, alldata)
}
#head(all_data)
#dim(all_data) # total_n = 18665 

nhanes_data <- all_data
save(nhanes_data, file="nhanes_data.RData")

