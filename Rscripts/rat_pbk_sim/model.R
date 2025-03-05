 # Rat Pyrethroid model
# Last modification: January 2019, Gina Song
# t is a vector of times output by the model; 
# State is the list of state variables in the model. The variables are as follows: ODOSE = oral dose,IVSWITCH = to put IV on or off, AORAL = Amount absorbed after oral, 
#                                                                                  AIV = amount absorbed via IV, AVE = amount in the vehicule compartment, 
#                                                                                  ALU = amount in the lumen, AFEC = amount in feces, ACYP = amount oxidized in liver by CYP
#                                                                                  ACAE = amount hydrolized by CES in liver, ACAEP = amount hydrolized by CES in plasma,
#                                                                                  AGI=amount in the GI,AL=amount in the liver, 
#                                                                                  AVF=amout in fat plasma, AF=amount in fat,ARP = amount in rapidly perfused tissues,
#                                                                                  AVSP = amount in slowly perfused tissues plasma, ASP = amount in slowly perfused tissues,
#                                                                                  AVB = amount in brain, AB = amount in blood,APLS = amount in plasma
# Parameters is the list of parameters for the model as describe in the scenarios files
# The output from the model are: "cliv" = conc in the liver, "cbrn"=conc in brain,"cpls"=conc in plasma, "odoseball" = dose mass balance, "bal"=mass balance,
#                               "vbal"=volumes mass balance,"qbal"=blood flow mass balance 


#source("GenericDiffusionLimited.R")
library(deSolve)
genericPyrethroidRatModel <- function(t,state,parameters) {
  with (
    
    as.list(c(state,parameters)),{
      
      ###############
      #concentrations
      ###############
      #Plasma arterial
      CPLS<-(APLS/VOLPLS)                       #plasma concentration (umol/L)
      CPLSf<-(APLS/VOLPLS)*FuPLS 		#unbound (free) concentration in plasma (umol/L)
      CPLSb<-(APLS/VOLPLS)*(1-FuPLS)	        #bound concentration in plasma  (umol/L)
      CPLSMETf<-(APLS/VOLPLS)*FuPLS/KMF   	#free concentration (umol/L) available for metabolism in plasma, adjusted for restricted clearance	
      CPLSng<-(APLS/VOLPLS)*MW  		#plasma concentration to be compared to the in vivo data (ng/ml)
      
      #GI
      CGI<-(AGI/VOLGI)				#GI concentration (umol/L)
      CGIf<-(AGI/VOLGI)*FuPLS/PGI		#concentration leaving GI venous plasma (umol/L)
      CVGI<-CPLSb + CGIf			#concentration in venous plasma leaving GI (bound+unbound)
            
      #Liver
      CL<-(AL/VOLLIVER)				#concentration in liver tissue (umol/L)
      CLf<-(AL/VOLLIVER)*FuPLS/PLIV		#free concentration leaving the liver tissue (umol/L)
      CLMETf<-(AL/VOLLIVER)*FuPLS/PLIV/KMF	#free concentration available for metabolism in liver adjusted for restricted clearance(umol/L)
      CVL<-CLf+CPLSb				#concentration in venous plasma leaving the liver (umol/L)
      CLIVng<-(AL/VOLLIVER)*MW  		#total liver concentration to be compared to the in vivo data (ng/ml)
      
      #Brain
      CVBf<-AVB/(VTBC*VOLBRAIN)			#concentration in brain plasma (umol/L)
      CBf<-(AB/((1-VTBC)*VOLBRAIN))*FuPLS/PBRN	#free concentration in brain tissue (umol/L), free fraction in brain=FuPLS/PBRN
      CVB<-CPLSb+CVBf				#concentration in venous plasma leaving brain (umol/L)
      CBRNng<-((AVB+AB)*MW)/VOLBRAIN		#total brain concentration to be compared to the in vivo data(ng/ml)
      
      #Fat
      CVFf<-AVF/(VTBC*VOLADIP)			#concentration in fat plasma (umol/L)
      CFf<-(AF/((1-VTBC)*VOLADIP))*FuPLS/PADIP  #free concentration in fat tissue (umol/L), free fraction in fat=FuPLS/PADIP
      CVF<-CPLSb+CVFf				#concentration in venous plasma leaving fat (umol/L)
      CFATng<-((AVF+AF)*MW)/VOLADIP		#total fat concentration to be compared to the in vivo data (ng/ml)
      
      #Rapidly perfused
      CRPf<-(ARP/VOLRP)*FuPLS/PRP		#free concentration leaving rapidly-perfused tissue (umol/L)
      CVRP<-CPLSb+CRPf				#concentration in venous plasma leaving rapidly-perfused tissue (umol/L)
      CRPng<-(ARP*MW)/VOLRP			#total fat concentration to be compared to the in vivo data (ng/ml)
      
      #Slowly-perfused
      CVSPf<-AVSP/(VTBC*VOLSP)			#concentration in slowly-perfused tissue plasma (umol/L)
      CSPf<-(ASP/((1-VTBC)*VOLSP))*FuPLS/PSP	#free concentration in slowly-perfused tissue (umol/L), free fraction in slowly-perfused tissue=FuPLS/PSP
      CVSP<-CVSPf+CPLSb				#concentration in venous plasma leaving slowly-perfused tissue (umol/L)
      CSPng<-((AVSP+ASP)*MW)/VOLSP		#total slowly-perfused tissue concentration to be compared to the in vivo data (ng/ml)
      
      #Plasma venous
      CV<-(QLIV*CVL+QADIP*CVF+QRP*CVRP+QSP*CVSP+QBRN*CVB)/CARDOUTP	#venous concentration (umol/L)
      
      ##########
      #Exposure
      ###########
      
      
      #Oral exposure
      RODOSE <- 0       # (umole)
      
      #IV exposure
      RIVSWTCH<-0       #switch for IV exposure (0 is off, 1 is on)
      
      RIV<-IVSWITCH*IDOSE/IVTLEN    #IV dose rate (umole/hour)
      
      ############
      #Equations
      ############
      
      #Vehicle compartment
      RVE<-  KLV*ALU-(KVL*AVE)-(KFEC*AVE)               #rate of change in vehicle (umol/hr)
      
      
      #GI lumen
      RLU<- KVL*AVE-KLV*ALU-KA*ALU-KFEC*ALU             #rate of change in lumen (umol/hr)
      
      #Fecal excretion 
      RFEC<- KFEC*(ALU+AVE)  				#rate of fecal excretion (umol/hr)
      
      #Absorption from lumen to the liver
      RORAL<- KA*ALU    				#rate of absorption (umol/hr)
      
	    #GI- Flow limited	
      RAGI<- QLIVGI*(CPLSf-CGIf)+(1-LYMPHSWTCH)*RORAL	#rate of change in GI tissue (umol/hr)
      
      #Liver- Flow limited
      #Liver metabolism
      RCYP<- (VCYP*CLMETf)/(KMCYP+CLMETf)  		                   #oxidation in liver by CYP (umol/hr)
      RCAE<- (VCAEM*CLMETf)/(KMCAEM+CLMETf)+(VCAEC*CLMETf)/(KMCAEC+CLMETf) #hydrolysis by CaE in liver (umol/hr)
      RAL<- QLIVH*CPLSf+QLIVGI*CGIf-QLIV*CLf-RCYP-RCAE	                   #rate of change in liver plasma (umol/hr)
      
      #Fat- Diffusion limited
      RAVF<- QADIP*(CPLSf-CVFf)+PAF*(CFf-CVFf)		        #rate of change in fat plasma (umol/hr)
      RAF<- PAF*(CVFf-CFf)					#rate of change in fat tissue (umol/hr)
      
      #Rapidly perfused-Flow limited
      RARP<- QRP*(CPLSf-CRPf)					#rate of change in rapidly-perfused tissue (umol/hr)
      
      #Slowly perfused-Diffusion limited
      RAVSP<- QSP*(CPLSf-CVSPf)+PASP*(CSPf-CVSPf)	#rate of change in slowly-perfused tissue plasma (umol/hr)
      RASP<- PASP*(CVSPf-CSPf)				#rate of change in slowly-perfused tissue (umol/hr)
      
      #Brain-Diffusion limited
      RAVB<- QBRN*(CPLSf-CVBf)+PAB*(CBf-CVBf)		#rate of change in brain plasma (umol/hr)
      RAB<- PAB*(CVBf-CBf)				#rate of change in brain tissue (umol/hr)
      
      #Plasma
      #Plasma metabolism
      RCAEP<- (VPLS*CPLSMETf)/(KMPLS+CPLSMETf) 		        #hydrolysis in plasma by CaE (umol/hr)
      RAPLS<- CARDOUTP*(CV-CPLS)-RCAEP+(LYMPHSWTCH*RORAL)+RIV	#rate of change in arterial plasma (umol/hr)
      
      #Calculation of clearance parameters
      CLINT_VIVO_ESTIMATED<-((VCYP/KMCYP)+(VCAEM/KMCAEM)+(VCAEC/KMCAEC))         #intrinsic clearance estimated directly based on in vitro measured Vmax and Km (L/hr) 
      CLINT_VIVO<-CLINT_VIVO_ESTIMATED/KMF   					 #intrinsic clearance in vivo adjusted for restrictive clearance (L/hr)
      CLH<- (CLINT_VIVO*QLIV)/((QLIV/FuPLS)+CLINT_VIVO)   			 #hepatic clearance (L/hr), assuming flow-limted, Yoon et al., TIV 2012, consistent with traditional Clh equation
      
      
      #Area Under Curves
      RAUCBRAIN<-0					#Area-under the concentration vs. time curve for brain(set to be 0)
      RAUCLIVER<-0					#Area-under the concentration vs. time curve for liver(set to be 0)
      RAUCADIPOSE<-0					#Area-under the concentration vs. time curve for fat(set to be 0)
      RAUCBLOOD<-0					#Area-under the concentration vs. time curve for blood(set to be 0)
      RAUCRBC<-0					#Area-under the concentration vs. time curve for RBC(set to be 0)
      
      #Mass balance
      TOTODOSE<-ODOSE-ALU-AVE-AORAL-AFEC  	 	         #mass balance for oral dose in vehicle, GI lumen and fecal excretion
      TOTDOSE<- AORAL + IDOSE 				         #total dose absorbed
      TOTBODY<-AVF+AF+ARP+AVSP+ASP+AVB+AB+APLS+AGI+AL	         #mass in body
      TOTCLEAR<-ACYP+ACAE+ACAEP				         #total amount cleared by metabolism in plasma and liver
      TMASS<-TOTDOSE-TOTBODY-TOTCLEAR			         #mass balance (should be 0)
      TVOL<-BODYWT*0.87 -(VOLLIVER+VOLADIP+VOLBLOOD+VOLSP+VOLRP+VOLBRAIN+VOLGI) #tissue volume balance (should be 0)
      TFLW<-CARDOUTP-(QADIP+QLIV+QBRN+QSP+QRP)		         # tissue plasma flow balance (should be 0)
      
      list(c(RODOSE, RIVSWTCH, RORAL, RIV, RVE, RLU, RFEC, RCYP, RCAE, RCAEP, RAGI, RAL, RAVF, RAF, RARP, RAVSP, RASP, RAVB, RAB, RAPLS, RAUCBRAIN, RAUCLIVER, RAUCADIPOSE, RAUCBLOOD, RAUCRBC),
           "cliv" = CLIVng,"cbrn"=CBRNng,"cpls"=CPLSng,"odosebal" = TOTODOSE,"bal"=TMASS,"vbal"=TVOL,"qbal"=TFLW) #,"cliv"=CLIV,"cl"=CL,"cvff"=CVFf,"cff"=CFf,"cfat"=CFAT,"crpf"=CRPF,"CRP"=CRP,
      
    }
  )
}
