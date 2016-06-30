#R Script creating a population of individuals with different covariates that Bayes and simulation methods will be tested on
#Source functions file
source(paste(work.dir,"Pop4_Inflix_Functions_File.R",sep = ""))
set.seed(123456)
#----------------------------------------------------------------------------------
#Define population's characteristics
#Only going to pre-specify weight as 70 kg and randomly generate ADA_TIME, BASE_ALB and FINAL_ALB
ID <- seq(from = 1,to = n,by = 1)	#Sequence of individual IDs
SIM <- sort(c(rep(seq(from = 0,to = nsim,by = 1),times = n)))	#Sequence of simulation identifiers
WT <- 70 #Weight, kg
AMP_ALB <-	0.1	#Amplitude for albumin sine wave
FREQ_ALB <- 1/60	#Frequency for albumin sine wave, number of oscillations per unit of time
PHASE_ALB <- 0	#Phase for albumin sine wave, where in its cycle the oscillation is a time = 0

#Create a data frame - one row per individual of covariate and random effects
ID.data <- data.frame(SIM,ID)
ID.data2 <- ID.data[ID.data$SIM != 0,]
cov.data <- ID.data

#Assign ID values to specific groups of covariate values
#ADA_TIME
ADA_TIME2	<- seq(from = 1,to = 3,by = 1)	#IDs with ADA present in the second sampling interval
ADA_TIME3	<- seq(from = 4,to = 6,by = 1)	#IDs with ADA present in the third sampling interval
ADA_TIME4	<- seq(from = 7,to = 9,by = 1)	#IDs with ADA present in the fourth sampling interval
ADA_TIME0	<- seq(from = 10,to = 12,by = 1)	#IDs with ADA never present

cov.data$ADA_TIME <- NA	#Add a ADA_TIME column
cov.data$ADA_TIME[cov.data$ID %in% ADA_TIME2] <- 154	#ADA in the second sampling interval
cov.data$ADA_TIME[cov.data$ID %in% ADA_TIME3] <- 294	#ADA in the third sampling interval
cov.data$ADA_TIME[cov.data$ID %in% ADA_TIME4] <- 462	#ADA in the fourth sampling interval
cov.data$ADA_TIME[cov.data$ID %in% ADA_TIME0] <- 646	#ADA never present (beyond the maximum of the TIME sequence)

#ALBUMIN
#BASE_ALB == 4 for all individuals
cov.data$BASE_ALB <- 4	#BASE_ALB == 4

#FINAL_ALB
FINAL_ALB3 <- c(1,4,7,10)	#ID's with FINAL_ALB == 3
FINAL_ALB4 <- c(2,5,8,11)	#ID's with FINAL_ALB == 4
FINAL_ALB5 <- c(3,6,9,12)	#ID's with FINAL_ALB == 5

cov.data$FINAL_ALB <- NA	#Add a FINAL_ALB column
cov.data$FINAL_ALB[cov.data$ID %in% FINAL_ALB3] <- 3	#FINAL_ALB == 3
cov.data$FINAL_ALB[cov.data$ID %in% FINAL_ALB4] <- 4	#FINAL_ALB == 4
cov.data$FINAL_ALB[cov.data$ID %in% FINAL_ALB5] <- 5	#FINAL_ALB == 5

#Simulate random effect parameters
#Simulating a baseline ETA and a final ETA to accommodate random changes in the individual that cannot be explained by model covariates
#Clearance stays as one value for ETA as it has a few time-dependent covariates on it
cov.data <- cov.data[with(cov.data, order(cov.data$ID,cov.data$SIM)), ]	#Sort by ID then SIM
cov.data$ETA1[cov.data$SIM != 0] <- rnorm(n*nsim,mean = 0,sd = PPVCL)	#ETA for clearance
cov.data$BASE_ETA2[cov.data$SIM != 0] <- rnorm(n*nsim,mean = 0,sd = PPVV1)	#Baseline ETA for V1
cov.data$FINAL_ETA2[cov.data$SIM != 0] <- log(exp(cov.data$BASE_ETA2[cov.data$SIM != 0])*c(0.9,1.1))	#Final ETA for V1
cov.data$BASE_ETA3[cov.data$SIM != 0] <- rnorm(n*nsim,mean = 0,sd = PPVQ)	#Baseline ETA for Q
cov.data$FINAL_ETA3[cov.data$SIM != 0] <- log(exp(cov.data$BASE_ETA3[cov.data$SIM != 0])*c(1.1,0.9))	#Final ETA for Q
cov.data$BASE_ETA4[cov.data$SIM != 0] <- rnorm(n*nsim,mean = 0,sd = PPVV2)	#Baseline ETA for V2
cov.data$FINAL_ETA4[cov.data$SIM != 0] <- log(exp(cov.data$BASE_ETA4[cov.data$SIM != 0])*c(0.9,1.1))	#Final ETA for V2

#For individuals in SIM == 0, make ETA == 0 (population predicted)
cov.data$ETA1[cov.data$SIM == 0] <- 0
cov.data$BASE_ETA2[cov.data$SIM == 0] <- 0
cov.data$FINAL_ETA2[cov.data$SIM == 0] <- 0
cov.data$BASE_ETA3[cov.data$SIM == 0] <- 0
cov.data$FINAL_ETA3[cov.data$SIM == 0] <- 0
cov.data$BASE_ETA4[cov.data$SIM == 0] <- 0
cov.data$FINAL_ETA4[cov.data$SIM == 0] <- 0

#----------------------------------------------------------------------------------
#Data frame of individual characteristics
pop.data <- lapply(cov.data,rep.int,times = length(TIME))
pop.data <- as.data.frame(pop.data)
pop.data <- pop.data[with(pop.data, order(pop.data$SIM,pop.data$ID)), ]	
pop.data$TIME <- TIME

#Calculate ETA values for all time-points
#A linear function containing the baseline ETA (BASE_ETA) and their final ETA (FINAL_ETA)
pop.data <- ddply(pop.data, .(SIM,ID), eta.function)

#Calculate albumin concentrations for each individual for all time-points
#A linear function containing the baseline albumin (BASE_ALB) and their last albumin (FINAL_ALB)
pop.data <- ddply(pop.data, .(SIM,ID), albumin.function)

#Flag if ADA are present for each individual for all time-points
#This assumes that once a person develops ADA, they stay with ADA
pop.data <- ddply(pop.data, .(SIM,ID), ada.function)

#Create a data frame of covariate values for each individual a key time-points
#i.e., baseline (0), day 98, 210, 378 and 546
cov.time.data <- subset(pop.data,select = c(ID,SIM,TIME,ALB,ADA))
cov.time.data <- cov.time.data[cov.time.data$TIME %in% c(0,98,210,378,546),]
 
#----------------------------------------------------------------------------------
##################
##_SIMULATION_1_##
##################
#In all scenarios, the initial compartment conditions are zero and the doses in the first sampling interval are 5 mg/kg
DOSE1 <- 5
pop.data$DOSE <- 0
pop.data$DOSE[pop.data$TIME %in% TIMEi] <- DOSE1
#Remove some unnecessary albumin, ADA_TIME and ETA columns, and reorder the others
pop.data <- pop.data[c("SIM","ID","TIME","DOSE","ADA","ALB","ETA1","ETA2","ETA3","ETA4")]
#Input data for the first simulation only contains time-points for the first sampling interval
input.data.sim1 <- pop.data[pop.data$TIME %in% TIME1all,]
#Assign compartment initial conditions in the input data frame
input.data.sim1$A1 <- 0	#Compartment 1 (central) initial conditions
input.data.sim1$A2 <- 0	#Compartment 2 (peripheral) initial conditions
input.data.sim1$AUT <- 0	#Time under the target trough concentration, 3 mg/L

#Apply simulate.conc.cmpf to each individual in input.data
sim.data1 <- ddply(input.data.sim1, .(SIM,ID), simulate.conc.cmpf)
sim.data1 <- as.data.frame(sim.data1)
sim.data1 <- sim.data1[sim.data1$TIME < max(TIME1all),]	#Remove the last time-point because some of my other code was lazy and this makes the plot look ugly
#This is just a time-point marking when the infusion administered just after sampling ended
