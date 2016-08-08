# Time-weighted Bayes project
# Script for simulating characteristics of the study population
# Does not simulation concentrations
# ------------------------------------------------------------------------------
# Source the other R scripts and execute
	work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/"
	source(paste0(work.dir,"functions.R"))

# ------------------------------------------------------------------------------
# Define population's characteristics
# Only going to pre-specify weight as 70 kg and randomly generate ADA_TIME, BASE_ALB and FINAL_ALB
	ID <- seq(from = 1,to = n,by = 1)	# Sequence of individual IDs
	SIM <- sort(c(rep(seq(from = 1,to = nsim,by = 1),times = n)))	# Sequence of simulation identifiers

	# First sine wave
		AMP_ALB1 <-	0.05	# Amplitude for albumin sine wave
		FREQ_ALB1 <- 1/100	# Frequency for albumin sine wave, number of oscillations per unit of time
		PHASE_ALB1 <- 0	# Phase for albumin sine wave, where in its cycle the oscillation is a time = 0
	# Second sine wave
		AMP_ALB2 <- 0.05
		FREQ_ALB2 <- 1/60
		PHASE_ALB2 <- 7
	# Third sine wave
		AMP_ALB3 <- 0.05
		FREQ_ALB3 <- 1/5
		PHASE_ALB3 <- 0

# Create a data frame - one row per individual of covariate and random effects
	ID.data <- data.frame(SIM,ID)
	cov.data <- ID.data	# Transfer to a data frame called "cov.data"

# Assign ID values to specific groups of covariate values
# WEIGHT
	cov.data$WT <- 70 # kg
# ADA_TIME
	ADA_TIME2	<- seq(from = 1,to = 5,by = 1)	# IDs with ADA present in the second sampling interval
	ADA_TIME3	<- seq(from = 6,to = 10,by = 1)	# IDs with ADA present in the third sampling interval
	ADA_TIME4	<- seq(from = 11,to = 15,by = 1)	# IDs with ADA present in the fourth sampling interval
	ADA_TIME0	<- seq(from = 16,to = 20,by = 1)	# IDs with ADA never present
	# Fill in ADA_TIME based on ID
		cov.data$ADA_TIME <- NA	# Add a ADA_TIME column
		cov.data$ADA_TIME[cov.data$ID %in% ADA_TIME2] <- 154	# ADA in the second sampling interval
		cov.data$ADA_TIME[cov.data$ID %in% ADA_TIME3] <- 294	# ADA in the third sampling interval
		cov.data$ADA_TIME[cov.data$ID %in% ADA_TIME4] <- 462	# ADA in the fourth sampling interval
		cov.data$ADA_TIME[cov.data$ID %in% ADA_TIME0] <- 646	# ADA never present (beyond the maximum of the TIME sequence)

# ALBUMIN
# BASE_ALB == 4 for all individuals
	cov.data$BASE_ALB <- 4	# BASE_ALB == 4
# FINAL_ALB
	FINAL_ALB1 <- c(1,6,11,16)	# ID's with FINAL_ALB == 1
	FINAL_ALB3 <- c(2,7,12,17)	# ID's with FINAL_ALB == 3
	FINAL_ALB4 <- c(3,8,13,18)	# ID's with FINAL_ALB == 4
	FINAL_ALB5 <- c(4,9,14,19)	# ID's with FINAL_ALB == 5
	FINAL_ALB7 <- c(5,10,15,20)	# ID's with FINAL_ALB == 7
	# Fill in final albumin values based on ID
		cov.data$FINAL_ALB <- NA	#Add a FINAL_ALB column
		cov.data$FINAL_ALB[cov.data$ID %in% FINAL_ALB1] <- 1	# FINAL_ALB == 1
		cov.data$FINAL_ALB[cov.data$ID %in% FINAL_ALB3] <- 3	# FINAL_ALB == 3
		cov.data$FINAL_ALB[cov.data$ID %in% FINAL_ALB4] <- 4	# FINAL_ALB == 4
		cov.data$FINAL_ALB[cov.data$ID %in% FINAL_ALB5] <- 5	# FINAL_ALB == 5
		cov.data$FINAL_ALB[cov.data$ID %in% FINAL_ALB7] <- 7	# FINAL_ALB == 7

# Simulate random effect parameters
# Simulating a baseline ETA and a final ETA to accommodate random changes in the individual that cannot be explained by model covariates
	cov.data <- cov.data[with(cov.data, order(cov.data$ID,cov.data$SIM)), ]	# Sort by ID then SIM
	# Clearance stays as one value for ETA as it has a few time-dependent covariates on it
		cov.data$ETA1 <- rnorm(n*nsim,mean = 0,sd = PPVCL)	# ETA for clearance
		cov.data$BASE_ETA2 <- rnorm(n*nsim,mean = 0,sd = PPVV1)	# Baseline ETA for V1
		cov.data$FINAL_ETA2 <- log(exp(cov.data$BASE_ETA2)*runif(n*nsim,min = 0.7,max = 1.3))	# Final ETA for V1
		cov.data$BASE_ETA3 <- rnorm(n*nsim,mean = 0,sd = PPVQ)	# Baseline ETA for Q
		cov.data$FINAL_ETA3 <- log(exp(cov.data$BASE_ETA3)*runif(n*nsim,min = 0.7,max = 1.3))	# Final ETA for Q
		cov.data$BASE_ETA4 <- rnorm(n*nsim,mean = 0,sd = PPVV2)	# Baseline ETA for V2
		cov.data$FINAL_ETA4 <- log(exp(cov.data$BASE_ETA4)*runif(n*nsim,min = 0.7,max = 1.3))	# Final ETA for V2

# ------------------------------------------------------------------------------
# Data frame of individual characteristics
	pop.data <- lapply(cov.data,rep.int,times = length(TIME))
	pop.data <- as.data.frame(pop.data)
	pop.data <- pop.data[with(pop.data, order(pop.data$SIM,pop.data$ID)), ]
	pop.data$TIME <- TIME

# Calculate albumin concentrations for each individual for all time-points
	pop.data <- ddply(pop.data, .(SIM,ID), albumin.function)

# Flag if ADA are present for each individual for all time-points
	pop.data <- ddply(pop.data, .(SIM,ID), ada.function)

# Calculate ETA values for all time-points
	pop.data <- ddply(pop.data, .(SIM,ID), eta.function)

# Simulating the same residual error
	pop.data$ERRPRO <- rnorm(length(pop.data$TIME),mean = 0,sd = ERRPRO)	# Proportional residual error

# ------------------------------------------------------------------------------
# Write pop.data to a .csv file
	pop.data.filename <- "population_characteristics.csv"
	write.csv(pop.data,file = pop.data.filename,na = ".",quote = F,row.names = F)
