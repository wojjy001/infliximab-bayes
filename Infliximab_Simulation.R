#Simulating Population
#R Script for simulating infliximab concentrations for a population
#------------------------------------------------------------------------------------------
#Remove all current objects in the workspace
rm(list=ls(all=TRUE))

#Load package libraries
library(deSolve)	#Differential equation solver
library(ggplot2)	#Plotting
library(plyr)	#Split and rearrange data, ddply function
library(grid)	#Plotting
library(compiler)	#Compile repeatedly-called functions

#------------------------------------------------------------------------------------------
#Set a directory for where plots can be saved (best where this R script is saved)
setwd("/Volumes/Prosecutor/PhD/InfliximabBayes/")

#Define a custom ggplot2 theme
theme_bw2 <- theme_set(theme_bw(base_size = 16)) 

#Function for calculating the median, and 2.5th and 97.5th percentiles for plotting simulation results
sumfuncx <- function(x) {
	stat1 <- median(x)
	stat2 <- quantile(x, probs=0.025, names=F) 
	stat3 <- quantile(x, probs=0.975, names=F)
	stat4 <- length(x)
	result <- c("median"=stat1, "low"=stat2, "hi"=stat3, "n"=stat4)
	result
}

#------------------------------------------------------------------------------------------
#Assign patient population characteristics
n <- 1	#Number of individuals to be simulated
ID <- seq(from = 1,to = n,by = 1)	#Sequence of individuals

#Patient characteristics
WT <- 70 #kg
ADA_TIME <- 300	#Onset of anti-drug antibodies, days
ALB <- 4	#Albumin

#Dosing and time information
TIMEi <- c(0,14,42,98,154,210,266,322)	#Infusion times (0, 2, 6 weeks and then every 8 weeks)
INFD <- round(2/24,digits = 3)	#Infusion duration, 2 hours (in days)
TIMEend <- TIMEi+INFD	#Mark when the infusion ends
TIMEr <- seq(from = 0,to = 365,by = 10)	#Make a time sequence at regular times
TIME <- unique(sort(c(TIMEi,TIMEend,TIMEr,365)))

AMT <- 5*WT #5 mg/kg as per AMH
RATE <- AMT/INFD	#Infusion rate, mg/day

#------------------------------------------------------------------------------------------
#Define the values for the model's population parameters
#THETAs
	#PK Parameters
	POPCL <- 0.294
	POPV1 <- 3.33
	POPQ <- 0.0719
	POPV2 <- 1.14
	#Covariate Effects
	WT_CL <- 0.614	#Effect of weight on clearance
	WT_V1 <- 0.691	#Effect of weight on V1
	WT_Q <- 1.1	#Effect of weight on Q
	WT_V2 <- 0.59	#Effect of weight on V2
	ALB_CL <- -1.17	#Effect of albumin on clearance
	ADA_CL <- 0.257	#Effect of anti-drug antibodies on clearance
	
#OMEGAs (SDs)
	#PK Parameters
	PPVCL <- 0.327
	PPVV1 <- 0.150
	PPVQ <- 1.10
	PPVV2 <- 0.799
	
#SIGMAs (SDs)
	#PK Parameters
	ERRPRO <- 0.419

#------------------------------------------------------------------------------------------
#Simulate random effect parameters
set.seed(123456)
ETA1 <- rnorm(n,mean = 0,sd = PPVCL)	#ETAs for CL
ETA2 <- rnorm(n,mean = 0,sd = PPVV1)	#ETAs for V1
ETA3 <- rnorm(n,mean = 0,sd = PPVQ)	#ETAs for Q
ETA4 <- rnorm(n,mean = 0,sd = PPVV2)	#ETAs for V2
	
#Make a data frame of individual parameter information to be entered into the differential equation solver - "input.data" - time-independent data
#Each individual needs one row in the data frame
input.data <- data.frame(ID,WT,ADA_TIME,ALB,ETA1,ETA2,ETA3,ETA4)
if (n > 1) {
	input.data$ETA1[input.data$ID == 1] <- 0
	input.data$ETA2[input.data$ID == 1] <- 0
	input.data$ETA3[input.data$ID == 1] <- 0
	input.data$ETA4[input.data$ID == 1] <- 0
}

#------------------------------------------------------------------------------------------
#Set up a ADA_TIME function
END <- max(TIME)+1	#Maximum value in the TIME sequence + 1
TIMEada <- c(0,ADA_TIME,END)
RATEada <- c(0,1,1)	
step.ada <- approxfun(TIMEada, RATEada, method="const")

#Create an "event.data" data frame (for albumin)
event.data <- data.frame(var = c(rep(3,times = length(TIME))),
												time = TIME,
												value = step.ada(TIME),
												method = "rep")
													
#------------------------------------------------------------------------------------------
#Set up infusion function
#Define the infusion (when it starts, when it finishes, and the rate) - this uses the "approxfun" function to make a "forcing function" for infusion rate in the differential equations
#The function needs to continue long after the infusion ends - specify an "end" time for the function
#Specify a vector that marks the infusion's time events
TIMEinf <- sort(c(TIMEi,TIMEi+INFD,END))	#Something happens at TIME = 0 and TIME = INFD
#Specify a vector marking the infusion's rates	
RATEinf <- c(rep(c(RATE,0),times = length(TIMEi)),0)	#At TIME = 0, RATE = RATE and TIME = INFD, RATE = 0 (i.e., infusion has ended, and continues to be zero until the end of the function) 
#Define an interpolation function that returns the rate when given a time - "const"
step.doseinf <- approxfun(TIMEinf, RATEinf, method="const")
																
#------------------------------------------------------------------------------------------
#Function containing differential equations for amount in each compartment
DES <- function(T, A, PAR) {
	#9 differential equations
	#2 x PK compartments, 1 x time-dependent variable(s)	
	dAdt <- vector(length = 3)
	
	################
	##_COVARIATES_##
	################
	#Time-independent covariates
	WT <- PAR[1]	#Weight, kg
	ALB <- PAR[2]	#Albumin
	#Time-dependent covariates
	ADA <- A[3]	#Presence of anti-drug antibodies, 0 = No, 1 = Yes
	dAdt[3] = 0 #Set rate to zero so covariate values don't change unless there is an event

	##############
	##_PK BLOCK_##
	##############
	#Define which values in the PAR vector (defined as PARAMETER.list below) are the ETA values for PK parameters
	ETA1 <- PAR[3]	#PPV for CL
	ETA2 <- PAR[4]	#PPV for V1
	ETA3 <- PAR[5]	#PPV for Q
	ETA4 <- PAR[6]	#PPV for V2
	
	#PK Parameters
	CL <- POPCL*((WT/70)^WT_CL)*((ALB/4)^ALB_CL)*(1+ADA_CL*ADA)*exp(ETA1)
	V1 <- POPV1*((WT/70)^WT_V1)*exp(ETA2)
	Q <- POPQ*((WT/70)^WT_Q)*exp(ETA3)
	V2 <- POPV2*((WT/70)^WT_V2)*exp(ETA4)
	
	#Apply the infusion function to the specified times
	RATEin <- step.doseinf(T)
	
	#Differential equations for PK
	dAdt[1] = RATEin -Q/V1*A[1] +Q/V2*A[2] -CL/V1*A[1]	#Central compartment
	dAdt[2] = Q/V1*A[1] -Q/V2*A[2]	#Peripheral compartment
		
	c(list(dAdt),"V1"=V1)
}

#Compile DES function - it's called by lsoda for each individual in the dataset	
DES.cmpf <- cmpfun(DES)

#------------------------------------------------------------------------------------------
#Function for simulating concentrations for the ith patient
simulate.conc <- function(input.data) {	
	#List of parameters from input for the differential equation solver			
	PARAMETER.list <- c(WT = input.data$WT,	#Weight
											ALB = input.data$ALB,	#Albumin
											ETA1 = input.data$ETA1,	#PPV for CL
											ETA2 = input.data$ETA2,	#PPV for V1
											ETA3 = input.data$ETA3,	#PPV for Q
											ETA4 = input.data$ETA4)	#PPV for V2
									
	#Set initial compartment conditions
	A1_0 <- 0	#Central compartment (PK)
	A2_0 <- 0	#Peripheral compartment (PK)
	ADA_0 <- event.data$value[event.data$var == 3 & event.data$time == 0]	#ADA
	
	#Vector it!
	A_0 <- c(A1 = A1_0,A2 = A2_0,ADA = ADA_0)
	
	#Run differential equation solver for simulated variability data	
	var.data <- ode(y = A_0,TIME,DES.cmpf,PARAMETER.list,events = list(data = event.data))
	var.data <- as.data.frame(var.data)	
}

#Compile simulate.conc function	- it's called by ddply for each individual in the dataset	
simulate.conc.cmpf <- cmpfun(simulate.conc)
#Apply simulate.conc.cmpf to each individual in par.data
#Maintain their individual value for V1 for later calculations
sim.data <- ddply(input.data, .(ID), simulate.conc.cmpf)
names(sim.data)[6] <- "V1"

#------------------------------------------------------------------------------------------
#Calculate individual predictions
sim.data$CONC <- sim.data$A1/sim.data$V1

#Add RUV
ERR <- rnorm(length(sim.data$CONC),mean = 0,sd = ERRPRO)
sim.data$DV <- sim.data$CONC*(1+ERR)

#------------------------------------------------------------------------------------------
#At each time-point ("time") in sim.data, calculate the median, and 5th and 95th percentiles for predictions
statsCONC <- ddply(sim.data, .(time), function(sim.data) sumfuncx(sim.data$CONC))
names(statsCONC)[c(2,3,4)] <- c("median","low","hi")

#Pick out some of the simulated observations to be "sampled"
TIMEobs <- 120
DVobs <- sim.data$DV[sim.data$time %in% TIMEobs]

#Plot CONC simulation results		
plotobj1 <- NULL
plotobj1 <- ggplot(statsCONC)
plotobj1 <- plotobj1 + geom_ribbon(aes(x = time,ymin = low,ymax = hi), alpha = 0.3)
plotobj1 <- plotobj1 + geom_line(aes(x = time,y = CONC),data = sim.data[sim.data$ID == 1,], colour="red")
plotobj1 <- plotobj1 + scale_y_continuous("Infliximab Concentration (mg/L)\n", breaks=c(0,20,40,60,80,100,120))
plotobj1 <- plotobj1 + scale_x_continuous("\nTime (days)")
print(plotobj1)

