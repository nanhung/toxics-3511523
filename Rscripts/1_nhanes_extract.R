# The function was modified by bayesmarker package
# The NHANES recently updated the data site (02/2025) 

library(readxl) # read_excel

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
