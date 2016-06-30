#R Script containing universal functions for simulating and fitting a population
ID <- NA	#Assign ID object
interval <- NA	#Assign sampling interval object
status <- NA	#Assign	simulation or Bayesian estimation status object
method <- NA	#Assign time-weighting method object
covariate <- NA	#Assign covariate method object
#----------------------------------------------------------------------------------
#Load package libraries
library(deSolve)	#Differential equation solver
library(ggplot2)	#Plotting
library(plyr)	#Split and rearrange data, ddply function
library(grid)	#Plotting
library(MASS)	#mvrnorm function
library(MBESS)	#cor2cov function
library(compiler)	#Compile repeatedly-called functions

trough.target <- 3	#Set the target trough concentration for dose optimisation
trough.upper <- 5	#Set upper bound for trough concentrations

#Custom ggplot2 theme
theme_bw2 <- theme_set(theme_bw(base_size = 14))

#Function for calculating the median, and 2.5th and 97.5th percentiles for plotting simulation results
CI95lo <- function(x) quantile(x,probs = 0.025)
CI95hi <- function(x) quantile(x,probs = 0.975)

#Function for taking the last row of each individual's profile
lastperID <- function(x) {tail(x,1)}

#----------------------------------------------------------------------------------
#Define time sequences
TIMEi <- c(0,14,42,98,154,210,266,322,378,434,490)	#Infusion times (0, 2, 6 weeks and then every 8 weeks) in days
INFD <- round(2/24,digits = 3)	#Infusion duration, 2 hours (in days)

TIME1 <- seq(from = 0,to = 98,by = 14)	#Specify the first sampling interval, days
TIME1i <- TIME1[TIME1 %in% TIMEi]	#Times in TIME1 that are infusion times
TIME1iend <- TIME1i+INFD	#Mark when the infusions end in the first sampling interval
TIME1all <- unique(sort(c(TIME1,TIME1i,TIME1iend)))
dose.time1 <- min(TIME1)	#Identify the first dosing time in first sampling interval
sample.time1 <- max(TIME1)	#Identify the first sampling time

TIME2 <- seq(from = 98,to = 210,by = 14) #Specify the second sampling interval, days
TIME2i <- TIME2[TIME2 %in% TIMEi]	#Times in TIME2 that are infusion times
TIME2iend <- TIME2i+INFD	#Mark when the infusions end in the second sampling interval
TIME2all <- unique(sort(c(TIME2,TIME2i,TIME2iend)))
dose.time2 <- min(TIME2)	#Identify the first dosing time in second sampling interval
sample.time2 <- max(TIME2)	#Identify the second sampling time

TIME3 <- seq(from = 210,to = 378,by = 14)	#Specify the third sampling interval
TIME3i <- TIME3[TIME3 %in% TIMEi]	#Times in TIME3 that are infusion times
TIME3iend <- TIME3i+INFD	#Mark when the infusions end in the third sampling interval
TIME3all <- unique(sort(c(TIME3,TIME3i,TIME3iend)))
dose.time3 <- min(TIME3)	#Identify the first dosing time in third sampling interval
sample.time3 <- max(TIME3)	#Identify the third sampling time

TIME4 <- seq(from = 378,to = 546,by = 14)	#Specify a fourth sampling interval
#This is just a time sequence for optimising the last dose - but this dose isn't actually given
TIME4i <- TIME4[TIME4 %in% TIMEi]	#Times in TIME4 that are infusion times
TIME4iend <- TIME4i+INFD	#Mark when the infusions end in the fourth sampling interval
TIME4all <- unique(sort(c(TIME4,TIME4i,TIME4iend)))
dose.time4 <- min(TIME4)	#Identify the first dosing time in fourth sampling interval
sample.time4 <- max(TIME4)	#Identify the fourth sampling time (though not actually sampled)

sample.times <- c(sample.time1,sample.time2,sample.time3,sample.time4)
TIME <- unique(sort(c(TIME1all,TIME2all,TIME3all,TIME4all)))	#Combine all TIME sequences together
END <- max(TIME)+100	#Make an object specifying a time-point long after the overall last time-point (for onset of ADA function)

#----------------------------------------------------------------------------------
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

#---------------------------------------------------------------------------------
#Function for calculating changes in random effects
#A linear function containing the baseline ETA (BASE_ETA) and their last ETA (FINAL_ETA)
eta.function <- function(input.data) {
	#ETA2
	TIMEeta2 <- c(min(input.data$TIME),max(input.data$TIME))
	RATEeta2 <- c(head(input.data$BASE_ETA2,1),head(input.data$FINAL_ETA2,1))
	step.eta2 <- approxfun(TIMEeta2,RATEeta2,method = "linear")	#Linear function
	input.data$ETA2 <- step.eta2(input.data$TIME)	#Apply function to every time-point
	#ETA3
	TIMEeta3 <- c(min(input.data$TIME),max(input.data$TIME))
	RATEeta3 <- c(head(input.data$BASE_ETA3,1),head(input.data$FINAL_ETA3,1))
	step.eta3 <- approxfun(TIMEeta3,RATEeta3,method = "linear")	#Linear function
	input.data$ETA3 <- step.eta3(input.data$TIME)	#Apply function to every time-point
	#ETA4
	TIMEeta4 <- c(min(input.data$TIME),max(input.data$TIME))
	RATEeta4 <- c(head(input.data$BASE_ETA4,1),head(input.data$FINAL_ETA4,1))
	step.eta4 <- approxfun(TIMEeta4,RATEeta4,method = "linear")	#Linear function
	input.data$ETA4 <- step.eta4(input.data$TIME)	#Apply function to every time-point
	as.data.frame(input.data)
}

#---------------------------------------------------------------------------------
#Function for calculating albumin concentrations for each individual for all time-points
#A linear function containing the baseline albumin (BASE_ALB) and their last albumin (FINAL_ALB)
albumin.function <- function(input.data) {
	TIMEalb <- c(min(input.data$TIME),max(input.data$TIME))
	RATEalb <- c(head(input.data$BASE_ALB,1),head(input.data$FINAL_ALB,1))
	step.alb <- approxfun(TIMEalb,RATEalb,method = "linear")	#Linear function
	input.data$ALB <- step.alb(input.data$TIME)*(1+AMP_ALB*sin(2*pi*FREQ_ALB*input.data$TIME+PHASE_ALB))	#Apply function to every time-point
	as.data.frame(input.data)
}

#Albumin.function for bayes fitting - does not assume sine wave - assume linear relationship between the two time-points
albumin.function.bayes <- function(input.data) {
	TIMEalb <- c(input.data$TIME[1],tail(input.data$TIME,1))
	RATEalb <- c(input.data$BASE_ALB[1],input.data$sample.ALB[1])
	step.alb <- approxfun(TIMEalb,RATEalb,method = "linear")	#Linearly extrapolate between the two observations
	input.data$ALB <- step.alb(input.data$TIME)	#Apply step.alb to each time-point in TIME
	input.data <- as.data.frame(input.data)
}

#----------------------------------------------------------------------------------
#Function for flagging if ADA are present for each individual for all time-points
#This assumes that once a person develops ADA, they stay with ADA
ada.function <- function(input.data) {
	TIMEada <- c(min(input.data$TIME),input.data$ADA_TIME[1],END)	#Specify times when ADA changes
	RATEada <- c(0,1,1)	#Specify the values for it to change to
	step.ada <- approxfun(TIMEada,RATEada,method = "const")	#Step function
	input.data$ADA <- step.ada(input.data$TIME)	#Apply function to every time-point
	as.data.frame(input.data)
}

#----------------------------------------------------------------------------------
#Function for collecting the cumulative time spent under the target trough concentration at each sampling time for each individual
AUT.collect.function <- function(input.data) {
	DOSE1 <- 5	#Time == 0
	AUT1 <- input.data$AUT[input.data$TIME == 14]	#AUT at the end of dosing interval
	DOSE2 <- 5	#Time == 14
	AUT2 <- input.data$AUT[input.data$TIME == 42]	#Cumulative AUT at the end of dosing interval
	DOSE3 <- 5	#TIME == 42
	AUT3 <- input.data$AUT[input.data$TIME == 98]	#Cumulative AUT at the end of dosing interval	
	DOSE4 <- input.data$DOSE[input.data$TIME == 98]	#TIME == 98
	AUT4 <- input.data$AUT[input.data$TIME == 154]	#Cumulative AUT at the end of dosing interval		
	DOSE5 <- input.data$DOSE[input.data$TIME == 154]	#TIME == 154
	AUT5 <- input.data$AUT[input.data$TIME == 210]	#Cumulative AUT at the end of dosing interval			
	DOSE6 <- input.data$DOSE[input.data$TIME == 210]	#TIME == 210
	AUT6 <- input.data$AUT[input.data$TIME == 266]	#Cumulative AUT at the end of dosing interval				
	DOSE7 <- input.data$DOSE[input.data$TIME == 266]	#TIME == 266
	AUT7 <- input.data$AUT[input.data$TIME == 322]	#Cumulative AUT at the end of dosing interval					
	DOSE8 <- input.data$DOSE[input.data$TIME == 322]	#TIME == 322
	AUT8 <- input.data$AUT[input.data$TIME == 378]	#Cumulative AUT at the end of dosing interval					
	DOSE9 <- input.data$DOSE[input.data$TIME == 378]	#TIME == 378
	AUT9 <- input.data$AUT[input.data$TIME == 434]	#Cumulative AUT at the end of dosing interval	
	DOSE10 <- input.data$DOSE[input.data$TIME == 434]	#TIME == 434
	AUT10 <- input.data$AUT[input.data$TIME == 490]	#Cumulative AUT at the end of dosing interval		
	DOSE11 <- input.data$DOSE[input.data$TIME == 490]	#TIME == 490	
	AUT11 <- input.data$AUT[input.data$TIME == 546]	#Cumulative AUT at the end of dosing interval		
	TOTAL_AMT <- sum(DOSE1,DOSE2,DOSE3,DOSE4,DOSE5,DOSE6,DOSE7,DOSE8,DOSE9,DOSE10,DOSE11)
	AUT.collect.data <- data.frame(ID = input.data$ID[1],
																SIM = input.data$SIM[1],
																DOSE1,AUT1,
																DOSE2,AUT2,
																DOSE3,AUT3,
																DOSE4,AUT4,
																DOSE5,AUT5,
																DOSE6,AUT6,
																DOSE7,AUT7,
																DOSE8,AUT8,
																DOSE9,AUT9,
																DOSE10,AUT10,
																DOSE11,AUT11,
																TOTAL_AMT)
}

#----------------------------------------------------------------------------------
#Function for simulating concentrations for the ith patient
simulate.conc <- function(input.data) {
	TIMEseq <- input.data$TIME #Pull out the specific time-sequence from the input.data
	#Create a "par.data" data frame (for non-time dependent parameters such at WT and ETAs)
	par.data <- head(input.data, 1)	#Collect individual information from the first time-point
	par.data <- par.data[c("ETA1","A1","A2","AUT")]

	#Set up "event.data" to specify changes in albumin, ADA status and dose
	event.data <- data.frame(var = c(
																	rep(3,times = length(TIMEseq)),	#ALB "compartment"
																	rep(4,times = length(TIMEseq)),	#ADA "compartment"
																	rep(5,times = length(TIMEseq)),	#ETA2 "compartment"
																	rep(6,times = length(TIMEseq)),	#ETA3 "compartment"
																	rep(7,times = length(TIMEseq)),	#ETA4 "compartment"
																	rep(8,times = length(TIMEseq))	#DOSE "compartment"
																	),
														time = c(TIMEseq,TIMEseq,TIMEseq,TIMEseq,TIMEseq,TIMEseq),
														value = c(input.data$ALB,
																			input.data$ADA,
																			input.data$ETA2,
																			input.data$ETA3,
																			input.data$ETA4,
																			input.data$DOSE),
														method = c(
																	rep("rep",times = length(TIMEseq)),
																	rep("rep",times = length(TIMEseq)),
																	rep("rep",times = length(TIMEseq)),
																	rep("rep",times = length(TIMEseq)),
																	rep("rep",times = length(TIMEseq)),
																	rep("rep",times = length(TIMEseq))
																	)
													)	#Brackets closing "data.frame"

	#List of parameters from input for the differential equation solver	(time-independent parameters)
	PARAMETER.list <- c(ETA1 = par.data$ETA1)	#ETA1

	#Function containing differential equations for amount in each compartment
	DES <- function(T, A, PAR) {
		#9 differential equations
		#2 x PK compartments, 6 x time-dependent variable(s), 1 x trough compartments
		dAdt <- vector(length = 9)

		################
		##_COVARIATES_##
		################
		#Time-dependent covariates
		ALB <- A[3] #Albumin
		ADA <- A[4]	#Anti-drug antibodies
		#Set rate to zero so covariate values don't change unless there is an event
		dAdt[3] = 0 #ALB
		dAdt[4] = 0	#ADA
		#Time-independent covariates
		WT <- 70	#Weight, kg

		##############
		##_PK BLOCK_##
		##############
		#Define which values in the PARAMETER.list vector are ETAs
		ETA1 <- PAR[1]	#PPV for CL
		ETA2 <- A[5]	#PPV for V1
		ETA3 <- A[6]	#PPV for Q
		ETA4 <- A[7]	#PPV for V2
		#Set rate to zero so ETA values don't change unless there is an event
		dAdt[5] = 0	#ETA2
		dAdt[6] = 0	#ETA3
		dAdt[7] = 0	#ETA4

		#Individual PK Parameters
		CL <- POPCL*((WT/70)^WT_CL)*((ALB/4)^ALB_CL)*(1+ADA_CL*ADA)*exp(ETA1)
		V1 <- POPV1*((WT/70)^WT_V1)*exp(ETA2)
		Q <- POPQ*((WT/70)^WT_Q)*exp(ETA3)
		V2 <- POPV2*((WT/70)^WT_V2)*exp(ETA4)

		#Dosing information
		DOSE <- A[8]	#Dose in mg/kg
		dAdt[8] = 0	#Rate of change in dose, doesn't change unless event is specified
		AMT <- DOSE*WT	#Amount in mg
		RATE <- AMT/INFD	#mg/day

		#Differential equations for PK
		dAdt[1] = RATE -Q/V1*A[1] +Q/V2*A[2] -CL/V1*A[1]	#Central compartment
		dAdt[2] = Q/V1*A[1] -Q/V2*A[2]	#Peripheral compartment
		
		############
		##_TROUGH_##
		############
		#Duration concentrations are below 3 mg/L
		dAdt[9] = 0
		CONC <- A[1]/V1 #Concentration in the central compartment
		if (CONC < 3) {dAdt[9] = 1}
		
		c(list(dAdt),"CL" = CL,"V1" = V1,"Q" = Q,"V2" = V2)
	}
	#Compile DES function - it's called by lsoda for each individual in the dataset
	DES.cmpf <- cmpfun(DES)

	#Set initial compartment conditions
	A1_0 <- par.data$A1	#Central compartment (PK)
	A2_0 <- par.data$A2	#Peripheral compartment (PK)
	ALB_0 <- event.data$value[event.data$var == 3 & event.data$time == head(event.data$time,1)]	#ALB
	ADA_0 <- event.data$value[event.data$var == 4 & event.data$time == head(event.data$time,1)]	#ADA
	ETA2_0 <- event.data$value[event.data$var == 5 & event.data$time == head(event.data$time,1)]	#ETA2
	ETA3_0 <- event.data$value[event.data$var == 6 & event.data$time == head(event.data$time,1)]	#ETA3
	ETA4_0 <- event.data$value[event.data$var == 7 & event.data$time == head(event.data$time,1)]	#ETA4
	DOSE_0 <- event.data$value[event.data$var == 8 & event.data$time == head(event.data$time,1)] #DOSE
	AUT_0 <- par.data$AUT	#Area under the target trough concentration (3 mg/L)
	A_0 <- c(A1 = A1_0,A2 = A2_0,ALB = ALB_0,ADA = ADA_0,ETA2 = ETA2_0,ETA3 = ETA3_0,ETA4 = ETA4_0,DOSE = DOSE_0,AUT = AUT_0)	#Vector of initial conditions

	#Run differential equation solver for simulated concentration data
	conc.data <- ode(y = A_0,TIMEseq,DES.cmpf,parms = PARAMETER.list,events = list(data = event.data),method = "lsoda")
	conc.data <- as.data.frame(conc.data)
	final.data <- data.frame(ID = input.data$ID,
													TIME = TIMEseq,
													DOSE = input.data$DOSE,
													AMT = input.data$DOSE*70,
													RATE = input.data$DOSE*70/INFD,
													WT = 70,
													ADA = input.data$ADA,
													ALB = input.data$ALB,
													ETA1 = input.data$ETA1,
													ETA2 = input.data$ETA2,
													ETA3 = input.data$ETA3,
													ETA4 = input.data$ETA4,
													A1 = conc.data$A1,
													A2 = conc.data$A2,
													CL = conc.data$CL.ALB,
													V1 = conc.data$V1.ETA2,
													Q = conc.data$Q.ETA3,
													V2 = conc.data$V2.ETA4,
													CONC = conc.data$A1/conc.data$V1.ETA2,
													AUT = conc.data$AUT)
	final.data$AMT[final.data$RATE == 0] <- 0
	final.data <- as.data.frame(final.data)
}

#Compile simulate.conc function	- it's called by ddply for each individual in the dataset
simulate.conc.cmpf <- cmpfun(simulate.conc)

#----------------------------------------------------------------------------------
#Create a data frame with the correct dose AND correct covariate values AND correct random effects (using pop.data)
#Simulate concentrations for the "TRUE" individual given the "optimised" dose
#However, if the Bayesian fitting said that the individual's trough concentration was in the target range, then continue with the previous dose
create.input.sim.optim.function <- function(dose.sampling,	#Dose at beginning of previous interval
																						sim.data.bayes,	#Simulate concentrations based on Bayesian fitted parameters
																						sample.time,	#Sample time of previous interval
																						optim.dose,	#Optimised dose based on parameters from the previous interval
																						final.sampling,	#Conditions and covariates at time of sampling in previous interval
																						TIMEperiod) {	#Time period to create new data frame on
	create.input.sim.optim.individual <- function(input.data) {
		#Covariate and random effect data
			cov.ind.data <- pop.data[pop.data$ID == input.data$ID & pop.data$SIM == input.data$SIM,]	#True random effect and covariate values
			cov.ind.data <- cov.ind.data[cov.ind.data$TIME %in% TIMEperiod,]	#Subset for just the second sampling interval
		#Dosing information
			OPTIM_DOSE1 <- optim.dose$DOSE_OPT1[optim.dose$ID == input.data$ID & optim.dose$SIM == input.data$SIM]	#Individual's first optimised dose
			OPTIM_DOSE2 <- optim.dose$DOSE_OPT2[optim.dose$ID == input.data$ID & optim.dose$SIM == input.data$SIM]
			OPTIM_DOSE3 <- optim.dose$DOSE_OPT3[optim.dose$ID == input.data$ID & optim.dose$SIM == input.data$SIM]
			PREV_DOSE <- dose.sampling$DOSE[dose.sampling$ID == input.data$ID & dose.sampling$SIM == input.data$SIM]	#Individual's previous dose
			DEC_CONC <- sim.data.bayes$CONC[sim.data.bayes$TIME == sample.time & sim.data.bayes$ID == input.data$ID & sim.data.bayes$SIM == input.data$SIM]#Individual Bayes estimated concentration at time of sampling, the "decision concentration"
			DOSE <- rep(PREV_DOSE,times = length(TIMEperiod))	#Define a "DOSE" variable as just the previous dose, unless it meets the conditions below
			TIMEinf <- TIMEperiod[TIMEperiod %in% TIMEi]	#TIMEs in TIMEperiod that are infusion times
			if (DEC_CONC < 3 | DEC_CONC >= 5 & TIMEperiod == sample.time) {
				DOSE[TIMEperiod == TIMEinf[1]] <- OPTIM_DOSE1
				DOSE[TIMEperiod == TIMEinf[2]] <- OPTIM_DOSE2
				DOSE[TIMEperiod == TIMEinf[3]] <- OPTIM_DOSE3
			}
		#Initial condition information
			A1 <- final.sampling$A1[final.sampling$ID == input.data$ID & final.sampling$SIM == input.data$SIM]	#A1
			A2 <- final.sampling$A2[final.sampling$ID == input.data$ID & final.sampling$SIM == input.data$SIM]	#A2
			AUT <- final.sampling$AUT[final.sampling$ID == input.data$ID & final.sampling$SIM == input.data$SIM]	#AUT
	
		#Final data frame
		result <- data.frame(cov.ind.data,A1,A2,AUT)
		result$DOSE <- DOSE
		result$DOSE[!c(result$TIME %in% TIMEi)] <- 0
		as.data.frame(result)
	}
	if (status == "simulated") {
		result <- ddply(ID.data, .(SIM,ID), create.input.sim.optim.individual)	#ID.data = input.data
	}
	if (status == "Bayesian estimated") {
		result <- ddply(ID.data2, .(SIM,ID), create.input.sim.optim.individual)	#ID.data = input.data	
	}
	result
}

#----------------------------------------------------------------------------------
#Function for optimising doses for an interval using maximum likelihood estimation
optimise.dose.function <- function(input.data) {
	# if (status == "simulated") {
		# print(paste("Optimising dose for Patient ",input.data$ID[1]," in Simulation ",input.data$SIM[1]," based on ",status," parameters for the ",interval," dosing interval.",sep = ""))
	# }
	# if (status == "Bayesian estimated") {
		# print(paste(method,covariate,": Optimising dose for Patient ",input.data$ID[1]," in Simulation ",input.data$SIM[1]," based on ",status," parameters for the ",interval," dosing interval.",sep = ""))
	# }
	TIMEinf <- input.data$TIME[input.data$TIME %in% TIMEi]	#TIMEs in the input dataset that are infusion times
	if (max(input.data$TIME) != 546) {
		TIMEinf <- TIMEinf[-NROW(TIMEinf)]
	}
	if (max(input.data$TIME) == 546) {
		TIMEinf <- TIMEinf
	}
	TIMEtrough <- TIMEinf+56	#TIMEs in the input dataset that a troughs
	initial.dose <- 10	#Initial parameter estimate, setting to standard dose of 5 mg/kg
	initial.err <- 0.001	#Initial estimated for "error"
	initial.par <- c(rep(initial.dose,times = length(TIMEinf)),initial.err)
	par <- initial.par
	optimise.interval.function <- function(par) {
		DOSE1fit <- par[1]	#First dose in sampling interval to be optimised
		DOSE2fit <- par[2]	#Second dose in sampling interval to be optimised
		input.data$DOSE[input.data$TIME == TIMEinf[1]] <- DOSE1fit	#Add into "input.data"
		input.data$DOSE[input.data$TIME == TIMEinf[2]] <- DOSE2fit		
		if (length(TIMEinf) <= 2) {
			#If there are only 2 doses in the sampling interval...
			ERRfit <- par[3]	#Parameter to be optimised, minimise the error
		}
		if (length(TIMEinf) > 2) {
			#If there are more than 2 doses in the sampling interval
			DOSE3fit <- par[3]	#Third dose in the sampling interval to be optimised
			input.data$DOSE[input.data$TIME == TIMEinf[3]] <- DOSE3fit
			ERRfit <- par[4]	#Parameter to be optimised, minimise the error
		}
		input.data$DOSE[!c(input.data$TIME %in% TIMEi)] <- 0	#When not an infusion time, make the DOSE = 0, so that it doesn't keep giving doses
		optim.dose.data <- ddply(input.data, .(), simulate.conc.cmpf)	#Simulate infliximab concentrations when given DOSEfit
		yhat <- optim.dose.data$CONC[optim.dose.data$TIME %in% TIMEtrough]
		postsigma <- yhat*ERRfit	#Posterior error
		res <- dnorm(trough.target,yhat,postsigma,log = T)	#Minimise the error between target trough (3 mg/L) and individual's actual trough
		objective <- -1*sum(res)
		objective
	}
	result <- optim(par,optimise.interval.function,hessian = FALSE,method = "L-BFGS-B",lower = c(rep(0.1,times = length(TIMEinf)),0.0001),upper = c(rep(500,times = length(TIMEinf)),1))	#Find the dose that minimise the sum of squared error#Find the dose that minimise the sum of squared error
	if (length(TIMEinf) <= 2) {
		result.data <- data.frame(DOSE_OPT1 = result$par[1],
															DOSE_OPT2 = result$par[2],
															DOSE_OPT3 = NA,
															ERR = result$par[3],	#Produce a data frame of the resultant dose and the minimised error
															OBJ = result$value[1])
	}
	if (length(TIMEinf) > 2) {
		result.data <- data.frame(DOSE_OPT1 = result$par[1],
															DOSE_OPT2 = result$par[2],
															DOSE_OPT3 = result$par[3],
															ERR = result$par[4],	#Produce a data frame of the resultant dose and the minimised error
															OBJ = result$value[1])	
	}
	result.data
}
optimise.dose.function.cmpf <- cmpfun(optimise.dose.function)

#----------------------------------------------------------------------------------
#Create a data frame fit for Bayesian estimation of a sampling interval
create.input.bayes.data.function <- function(covariate,dose.sampling,init.sampling,final.sampling,dose.time,sample.time,bayes.result,TIME,method) {
	create.input.bayes.data.individual <- function(input.data) {
		DOSE <- dose.sampling$DOSE[dose.sampling$ID == input.data$ID & dose.sampling$SIM == input.data$SIM]	#Source dose from previous sampling interval
		A1 <- 0	#Initial conditions are zero - will be estimating right from the beginning since first dose
		A2 <- 0
		AUT <- 0
		DV <- final.sampling$DV[final.sampling$ID == input.data$ID & final.sampling$SIM == input.data$SIM]	#Source sample from previous sampling interval
		ADA <- cov.time.data$ADA[cov.time.data$TIME == sample.time & cov.time.data$ID == input.data$ID & cov.time.data$SIM == input.data$SIM]	#Source ADA status at the end of previous sampling interval
		BASE_ALB <- cov.time.data$ALB[cov.time.data$TIME == dose.time & cov.time.data$ID == input.data$ID & cov.time.data$SIM == input.data$SIM]	#Source baseline albumin from previous sampling interval
		FINAL_ALB <- cov.time.data$ALB[cov.time.data$TIME == sample.time & cov.time.data$ID == input.data$ID & cov.time.data$SIM == input.data$SIM]	#Source final albumin from previous sampling interval
		TIMEalb <- c(dose.time,sample.time)
		RATEalb <- c(BASE_ALB,FINAL_ALB)
		step.alb <- approxfun(TIMEalb,RATEalb,method = "linear")	#Linearly extrapolate between the two time-points
		ALB <- step.alb(TIME)		#Linearly extrapolate between the two time-points
		#Below depends on the "covariate" variable
		if (covariate == "NoADA") {	#Make ADA the population typical, i.e., 0
			ADA <- 0
		}
		if (covariate == "NoALB") {	#Make ALB the population typical, i.e., 4
			ALB <- 4
		}
		if (covariate == "NoCov") {	#Make both ADA and ALB the population typical
			ADA <- 0
			ALB <- 4
		}		
		ETA1 <- bayes.result$ETA1[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]	#Source estimated ETA value from the previous interval - this will serve as our initial estimate
		ETA2 <- bayes.result$ETA2[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]	
		ETA3 <- bayes.result$ETA3[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]	
		ETA4 <- bayes.result$ETA4[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]	
		result <- data.frame(TIME,DOSE,ETA1,ETA2,ETA3,ETA4,ADA,ALB,A1,A2,AUT,DV,method)
		result$DOSE[!c(result$TIME %in% TIMEi)] <- 0 #Reset non-infusion time dose values
		result$DV[result$TIME != sample.time] <- NA	#Remove DV value from non-sampling times
		result
	}
	result <- ddply(ID.data2, .(SIM,ID), create.input.bayes.data.individual)	#ID.data = input.data
}

#----------------------------------------------------------------------------------
#Create a data frame for simulating a concentration time-profile for each individual based on fitted individual parameters and "known" covariate values
create.input.bayesfit.data.function <- function(input.bayes.data,bayes.result) {
	create.input.bayesfit.data.individual <- function(input.data) {
		result <- input.bayes.data[input.bayes.data$ID == input.data$ID & input.bayes.data$SIM == input.data$SIM,]
		result <- result[-ncol(result)]	#Remove the DV column
		ETA1 <- bayes.result$ETA1[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]	#Add Bayes fitted ETA1
		ETA2 <- bayes.result$ETA2[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]
		ETA3 <- bayes.result$ETA3[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]
		ETA4 <- bayes.result$ETA4[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]
		result <- data.frame(result,ETA1,ETA2,ETA3,ETA4)
	}
	result <- ddply(ID.data2, .(SIM,ID), create.input.bayesfit.data.individual)	#ID.data = input.data
}

#----------------------------------------------------------------------------------
#Create a data frame for dose optimisation based on "known" covariate values, Bayesian estimated fitted parameters and start from conditions from the Bayesian simulation
create.input.optim.data.function <- function(bayesian.sampling,bayes.result,TIME) {
	create.input.optim.data.individual <- function(input.data) {
		ALB <- bayesian.sampling$ALB[bayesian.sampling$ID == input.data$ID & bayesian.sampling$SIM == input.data$SIM]
		ADA <- bayesian.sampling$ADA[bayesian.sampling$ID == input.data$ID & bayesian.sampling$SIM == input.data$SIM]
		A1 <- bayesian.sampling$A1[bayesian.sampling$ID == input.data$ID & bayesian.sampling$SIM == input.data$SIM]
		A2 <- bayesian.sampling$A2[bayesian.sampling$ID == input.data$ID & bayesian.sampling$SIM == input.data$SIM]
		CONC <- bayesian.sampling$CONC[bayesian.sampling$ID == input.data$ID & bayesian.sampling$SIM == input.data$SIM]
		AUT <- bayesian.sampling$AUT[bayesian.sampling$ID == input.data$ID & bayesian.sampling$SIM == input.data$SIM]
		ETA1 <- bayes.result$ETA1[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]	#Add Bayes fitted ETA1
		ETA2 <- bayes.result$ETA2[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]
		ETA3 <- bayes.result$ETA3[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]
		ETA4 <- bayes.result$ETA4[bayes.result$ID == input.data$ID & bayes.result$SIM == input.data$SIM]
		result <- data.frame(TIME,ALB,ADA,A1,A2,CONC,AUT,ETA1,ETA2,ETA3,ETA4)
	}
	if (status == "simulated") {
		result <- ddply(ID.data, .(SIM,ID), create.input.optim.data.individual)	#ID.data = input.data
	}
	if (status == "Bayesian estimated") {
		result <- ddply(ID.data2, .(SIM,ID), create.input.optim.data.individual)	#ID.data = input.data
	}
	result
}

#----------------------------------------------------------------------------------
#Function for estimating individual parameters (ETAs)
bayesian.function <- function(input.data) {
	# print(paste(method,covariate,": Bayesian estimating parameters for Patient ",input.data$ID[1]," in Simulation ",input.data$SIM[1]," for the ",interval," dosing interval.",sep = ""))
	#Initial parameter estimates
	ETA1init <- input.data$ETA1[1]
	ETA2init <- input.data$ETA2[1]
	ETA3init <- input.data$ETA3[1]
	ETA4init <- input.data$ETA4[1]
	initial.par <- c(exp(ETA1init),exp(ETA2init),exp(ETA3init),exp(ETA4init))	#Population values
	par <- initial.par
	#Observation - posterior
	Yobs <- input.data$DV	#Most of this will be NA except for the samples
		#Function for minimising the Bayes Objective Function Value
		bayes.ofv <- function(par) {
			#Assign initial parameter values to an ETA to be estimated
			ETA1fit <- log(par[1])	#Fit the exp(ETA) for CL
			ETA2fit <- log(par[2])	#Fit the exp(ETA) for V1
			ETA3fit <- log(par[3])	#Fit the exp(ETA) for Q
			ETA4fit <- log(par[4])	#Fit the exp(ETA) for V2
			#Make a data frame of individual-specific parameters required as input for the differential equation solver
			#Needs to be in the same format that was used when simulating the patient
			input.bayes.data <- input.data[c("SIM","ID","TIME","DOSE","ADA","ALB","A1","A2","AUT")]
			input.bayes.data$ETA1 <- ETA1fit
			input.bayes.data$ETA2 <- ETA2fit
			input.bayes.data$ETA3 <- ETA3fit
			input.bayes.data$ETA4 <- ETA4fit
			#Run the differential equation solver to obtain concentrations
			bayes.data <- ddply(input.bayes.data, .(), simulate.conc.cmpf)
			#Pull out the individual predicted concentrations
			Yhat <- bayes.data$CONC
			#If Yobsx was NA, then Yhat needs to be NA too (for calculating the log-likelihood)
			Yhat[is.na(Yobs) == T] <- NA
			#Posterior component (from the data)
			#Log densities of residuals
			#Residual error model, Y = IPRE*(1+ERR), Y = IPRE + IPRE*ERR
			TIMET <- max(bayes.data$TIME) - bayes.data$TIME	#Time since last observation
			method <- input.data$method[1]
			if (method == "NTimeWeight") {	#No time-weighting
				loglikpost.sd <- na.omit(Yhat)*ERRPRO
			}
			if (method == "Peck1.005") {	#Peck method, Q = 1.005
				loglikpost.sd <- na.omit(Yhat)*ERRPRO*1.005^TIMET[is.na(Yhat) == F]
			}
			if (method == "Peck1.01") {	#Peck method, Q = 1.01
				loglikpost.sd <- na.omit(Yhat)*ERRPRO*1.01^TIMET[is.na(Yhat) == F]
			}
			loglikpost <- dnorm(na.omit(Yobs),mean = na.omit(Yhat),sd = loglikpost.sd,log = T)
			#Prior component (from the model)
			ETA <- c(ETA1fit,ETA2fit,ETA3fit,ETA4fit)
			ETABSV <- c(PPVCL,PPVV1,PPVQ,PPVV2)
			loglikprior <- dnorm(ETA,mean = 0,sd = ETABSV,log = T)
			#Calculate the combined likelihood
			OFVBayes <- -1*sum(loglikpost,loglikprior)
			OFVBayes
		}
	#Optimise the ETA parameters to minimise the OFVBayes
	resultfit <- optim(par,bayes.ofv,hessian = FALSE,method = "L-BFGS-B",lower = c(0.001,0.001,0.001,0.001),upper = c(Inf,Inf,Inf,Inf),control = list(parscale = par,factr = 1e7))
	#Put results in a data frame - not essentially necessary, but if I use this for a population, this is the code that I will need to bind everyone's results together
	resultfit.data <- data.frame(ETA1 = log(resultfit$par[1]),
															ETA2 = log(resultfit$par[2]),
															ETA3 = log(resultfit$par[3]),
															ETA4 = log(resultfit$par[4]))
	resultfit.data
}
bayesian.function.cmpf <- cmpfun(bayesian.function)
