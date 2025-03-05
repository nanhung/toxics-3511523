# This code is to simulate Cis-PM dose matrics in PND90. 
# Code update from Mallick et al. (2020)

#clear workspace
rm(list = ls())

# Package
library(deSolve)

# Model and paramFile
source("Rscripts/rat_pbk_sim/model.R")
source("Rscripts/rat_pbk_sim/CPM_params_Oral_90d.R")


# ORAL and IV
params[["PDOSE"]] <-39.17 * 0.4 # BBMD BMDL 
params[["IVDOSE"]]<-0   # mg/kg

#set simulation Parameters
params[["TSTART"]] <-0    #start time for simulation (day)
params[["TOTDAYS"]] <-1   #end of simulation (day)

initial_params <- within(as.list(params),{
  # PARAMETER SCALING
  # TISSUE VOLUME CHECK
  # fractional slowly-perfused tissue volume
  VOLSPC<- 0.87-(VOLBLOODC + VOLADIPC+ VOLBRAINC + VOLLIVERC +VOLRPC+VOLGIC)	
  # total fractions of tissue volume (should be 1)
  VOLTOTALC <- VOLBLOODC + VOLADIPC+ VOLBRAINC + VOLLIVERC +VOLRPC+VOLSPC+VOLGIC	
  # SCALED TISSUE VOLUMES (L)
  VOLBLOOD <- VOLBLOODC * BODYWT    # L; BLOOD
  VOLADIP <- VOLADIPC * BODYWT      # L; ADIPOSE TISSUE
  VOLBRAIN <- VOLBRAINC * BODYWT    # L; BRAIN TISSUE
  VOLLIVER <- VOLLIVERC *BODYWT     # L; LIVER TISSUE
  VOLRP <-  VOLRPC*BODYWT           #L ; RICHLY PERFUSED
  VOLSP <- VOLSPC*BODYWT            #L ; SLOWLY PERFUSED
  VOLGI <- VOLGIC*BODYWT            #L ; GI
  VOLPLS <- VOLBLOOD*(1-HCT)        #L ; PLASMA
  VOLRBC <- VOLBLOOD*HCT            #L ; RED BLOOD CELL
  #PLASMA FLOW CHECK
  #fractional rapidly-perfused tissue plasma flow
  FRRPC<-1-(FRADIPC + FRBRNC  + FRLIVC +FRSPC)
  #total fractions of tissue plasma flow (should be 1)
  FRTOTALC <- FRADIPC + FRBRNC  + FRLIVC +FRRPC+FRSPC	
  # SCALED PLASMA FLOWS (L / HR)
  CARDOUTP <- CARDOUTPC            #L/HR; CARDIAC PLASMA OUTPUT 
  # Note that rat QC is plasma CO, No need HCT
  QADIP <- FRADIPC * CARDOUTP      #L/HR; ADIPOSE PLASMA FLOW
  QBRN <- FRBRNC *CARDOUTP         #L/HR; BRAIN PLASMA FLOW
  QLIV <- FRLIVC * CARDOUTP        #L/HR; LIVER TOTAL PLASMA FLOW
  QLIVH <- QLIV*FRLIVH             #L/HR; FRACTION OF LIVER ARTERIAL PLASMA FLOW
  QLIVGI <- QLIV*(1-FRLIVH)        #L/HR; FRACTION OF LIVER GI PLASMA FLOW
  QRP <- FRRPC* CARDOUTP           #L/HR; FRACTION OF RAPIDLY-PERFUSED TISSUE PLASMA FLOW
  QSP <- FRSPC*CARDOUTP            #L/HR; FRACTION OF SLOWLY-PERFUSED TISSUE PLASMA FLOW
  #Permiability constants
  PAF <- PAFC* (VOLADIP ^0.75)	  # fat (L/hr/Kg^0.75)
  PASP <- PASPC* (VOLSP ^0.75)        # Slowly-perfused tissue (L/hr/Kg^0.75)
  PAB <- PABC* (VOLBRAIN ^0.75)       # brain (L/hr/Kg^0.75)
  #Pyrethroid metabolism
  VCYP <-VCYPC  * VOLLIVER   #Vmax liver CYP scaled to whole liver (umol/hr)
  VCAEM<-VCAEMC * VOLLIVER   #Vmax liver microsomal-CES scaled to whole liver(umol/hr)
  VCAEC<-VCAECC * VOLLIVER   #Vmax liver cytosolic-CES scaled to whole liver (umol/hr)
  VPLS<-VPLSC * VOLPLS       #Vmax Plasma CES scaled to plasma(umol/hr)
  TSTOP <- TOTDAYS*24     #End of simulation (hours)
})

#calculate in vivo hepatic clearance
VCYP<-initial_params[["VCYP"]]			#Vmax liver CYP (umol/hr)
VCAEM<-initial_params[["VCAEM"]]		#Vmax liver microsomal-CES (umol/hr)
VCAEC<-initial_params[["VCAEC"]]		#Vmax liver cytosolic-CES (umol/hr)
KMCYP<-initial_params[["KMCYP"]]		#Km Liver CYP (umol/L)
KMCAEM<-initial_params[["KMCAEM"]]		#Km liver microsomal-CES (umol/L)
KMCAEC<-initial_params[["KMCAEC"]]		#Km liver cytosolic-CES (umol/L)
KMF  <- initial_params[["KMF"]]			#empirical adjustment factor for restricted clearance
QLIV <- initial_params[["QLIV"]]		#plasma flow to liver (L)
FuPLS <- initial_params[["FuPLS"]]		#unbound fraction

CLINT_VIVO_ESTIMATED<-((VCYP/KMCYP)+(VCAEM/KMCAEM)+(VCAEC/KMCAEC)) # (L/h) estimated in vivo Clint based on measured in vitro Clint
CLINT_VIVO<-CLINT_VIVO_ESTIMATED/KMF   # (L/h)  in vivo Clint adjusted for restrictive CL
CLH<- (CLINT_VIVO*QLIV)/((QLIV/FuPLS)+CLINT_VIVO) # (L/h) hepatic clearance


#initialize state variables
state <- c(ODOSE=0, IVSWITCH=0, AORAL=0, AIV=0, AVE=0, ALU=0, AFEC=0, 
  ACYP=0, ACAE=0, ACAEP=0, AGI=0, AL=0, AVF=0, AF=0, ARP=0, AVSP=0, ASP=0,
  AVB=0, AB=0, APLS=0, AUCBRAIN=0, AUCLIVER=0, AUCADIPOSE=0, AUCBLOOD=0, AUCRBC=0)
tstop = initial_params[["TSTOP"]]

#Define Events 
bw <- initial_params[["BODYWT"]]	#Body weight (kg)
mw <- initial_params[["MW"]]		#molecular weight (g/mole)
pdose <- initial_params[["PDOSE"]]	#oral dose (mg/kg)
ivdose <- initial_params[["IVDOSE"]]	#IV dose (mg/kg)
ivtlen <- initial_params[["IVTLEN"]]	#IV duration (hr)

# calculate IDOSE to use for mass balance.
initial_params[["IDOSE"]]<-(ivdose*1000*bw/mw)     #IV dose (umole)

if (pdose > 0){
  #var to change
  state_var = c("AVE","ODOSE")		#Mass in vehicle (umole)
  #Value  of change
  change_val1= (pdose*bw*1000/mw)		#oral dose (umole)
  change_val2= change_val1
  change_arr = c(change_val1,change_val2)
  #operation of event
  operation = c("rep","rep") 
  event_times <- c(0,0) 
  eventDat <- data.frame( 
    var = state_var,# rep(x = state_Var,each = length(event_times)),
    time = event_times,#rep(event_times,length(state_Var)),
    value = change_arr,#rep(x = change_arr,each = length(event_times)),
    method = operation#rep(x = operation,each = length(event_times))   
  )
}else if (ivdose >0){
  # var to change
  state_var1 = "IVSWITCH"
  state_var2 ="IVSWITCH"
  # Value  of change 
  change_val1= 1
  change_val2= 0
  # operation of event
  operation1 = "rep"
  operation2 = "rep"
  # times of event
  time1 = 0
  time2 = ivtlen
  eventDat <- data.frame(
    var = c(state_var1,state_var2),
    time = c(time1,time2),
    value = c(change_val1,change_val2),
    method = c(operation1,operation2)
  )
}
t_start <- as.numeric(0)
t_dur <- as.numeric(tstop)
times <- seq(t_start,t_dur,by=0.01)

#run the model
modelOutput<- ode(y = state, times = times,method = "lsodes",
                  func = genericPyrethroidRatModel, parms = initial_params,
                  events = list(data = eventDat))
result <- as.data.frame(modelOutput)

# Systematic POD
auc_plasma <- sum(diff(result[, "time"]) * (result[-1, "cpls"] + result[-length(result), "cpls"]) / 2) 
auc_brain <- sum(diff(result[, "time"]) * (result[-1, "cbrn"] + result[-length(result), "cbrn"]) / 2) 
dosemetrics <- data.frame(organ=c('Plasma', 'Brain'),
  Cmax=c(max(result$cpls), max(result$cbrn)),
  AUC=c(auc_plasma, auc_brain))
dosemetrics
