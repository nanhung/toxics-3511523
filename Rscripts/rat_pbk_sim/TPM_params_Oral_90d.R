# File created by Gina Song
# Last update by Marjory Moreau on March 2019
# Description of the parameters used in the PBPK model for TPM for PND 90 rats.

params <- c(
  BODYWT =0.3517, 	#body weight (kg)
  MW = 391.3,       	#molecular weight for DLM  (g/mole)  
  FuPLS=0.2,      	#unbound fraction in plasma
  LYMPHSWTCH =1,  	#switch for lymphatic absorption 
  KMF = 10,       	#empirical adjustment factor for decreased free concentration
  KA = 0.31,      	#absorption rate constant from lumen (1/hr)
  KLV=0.000026,		#Rate constant for pyrethroid transfer from gut lumen compartment to corn oil compartment (1/hr)
  KVL=0.20,		#Rate constant for pyrethroid transfer from corn oil compartment to gut lumen compartment (1/hr)
  KFEC=0.025,     	#Fecal excretion rate constant (1/hr)
  PDOSE =0.0,     	#oral dose (mg/kg)
  IVDOSE = 0.0,   	#IV dose (mg/kg)
  IVTLEN=0.01,    	#IV duration (hr)
  VOLBLOODC =  0.074,   #fraction of bodyweight to blood volume 
  VOLADIPC =  0.0513,   #fraction of bodyweight to fat 
  VOLBRAINC = 0.005587, #fraction of bodyweight to brain
  VOLLIVERC = 0.0387,   #fraction of bodyweight to liver
  VOLRPC = 0.0155,  	#fraction of bodyweight to rapidly-perfused tissue
  VOLGIC = 0.0437,  	#fraction of bodyweight to GI
  VTBC = 0.05,      	#fraction of blood volume in diffusion-limited tissues
  CARDOUTPC =3.1197,    #total plasma flow (L/hr)
  HCT = 0.45,     	#hematocrit
  FRADIPC = 0.07, 	#fraction of total plasma flow to fat
  FRBRNC = 0.02,  	#fraction of total plasma flow to brain
  FRLIVC = 0.183, 	#fraction of total plasma flow to liver (portal + hepatic arterial blood)
  FRLIVH = 0.05,  	#fraction of total plasma flow to hepatic arterial 
  FRSPC = 0.15,   	#fraction of total plasma flow to slowly-perfused tissue
  PBRN= 0.44,     	#brain PC
  PADIP= 68.7,    	#fat PC
  PLIV = 0.6,    	#liver PC
  PRP = 0.6,     	#rapidly-perfused tissue PC
  PSP  =3.94,     	#slowly-perfused tissue PC
  PGI = 0.6,     	#GI PC
  PAFC=1.5,      	#fat permeability-area cross product (L/hr/kg^0.75)
  PASPC=0.05,    	#slowly-perfused tissue permeability-area cross product (L/hr/kg^0.75)
  PABC =0.095,    	#brain permeability-area cross product (L/hr/kg^0.75)
  VCYPC = 26049.6, 	#Vmax for liver CYP  (umol/hr/kg)
  VCAEMC=28823.0,    	#Vmax for liver microsomal CES (umol/hr/kg)
  VCAECC=6902.5,   	#Vmax for liver cytosolic CES (umol/hr/kg)
  VPLSC=2667.6,    	#Vmax for plasma CES (umol/hr/kg)
  KMCYP=7.2,     	#Km for liver CYP (umol/l)
  KMCAEM=2.78,    	#Km for liver microsomal CES (umol/l)
  KMCAEC=0.42,     	#Km for liver cytosolic CES (umol/l)
  KMPLS=0.6     	#Km for plasma CES (umol/l)
)   
