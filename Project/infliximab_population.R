# Time-weighted Bayes project
# Script for simulating characteristics of the study population
# Does not simulation concentrations
#-------------------------------------------------------------------------------
# Source the other R scripts and execute
	# infliximab_functions.R
	source(paste0(work.dir,"infliximab_functions.R"))

# ------------------------------------------------------------------------------
# Define population's characteristics
# Only going to pre-specify weight as 70 kg and randomly generate ADA_TIME, BASE_ALB and FINAL_ALB
	n <- 12	# Number of seed individuals (where each seed individual has a different set of covariate values)
	nsim <- 1	# Number of simulations of the seed individuals to perform
	ID <- seq(from = 1,to = n,by = 1)	# Sequence of individual IDs
	SIM <- sort(c(rep(seq(from = 0,to = nsim,by = 1),times = n)))	# Sequence of simulation identifiers
	WT <- 70 # Weight, kg
	AMP_ALB <-	0.1	# Amplitude for albumin sine wave
	FREQ_ALB <- 1/60	# Frequency for albumin sine wave, number of oscillations per unit of time
	PHASE_ALB <- 0	# Phase for albumin sine wave, where in its cycle the oscillation is a time = 0

# Create a data frame - one row per individual of covariate and random effects
	ID.data <- data.frame(SIM,ID)
	ID.data2 <- ID.data[ID.data$SIM != 0,]	# SIM = 0 is reserved for population typical
	cov.data <- ID.data	# Transfer to a data frame called "cov.data"

# Assign ID values to specific groups of covariate values
# ADA_TIME
	ADA_TIME2	<- seq(from = 1,to = 3,by = 1)	# IDs with ADA present in the second sampling interval
	ADA_TIME3	<- seq(from = 4,to = 6,by = 1)	# IDs with ADA present in the third sampling interval
	ADA_TIME4	<- seq(from = 7,to = 9,by = 1)	# IDs with ADA present in the fourth sampling interval
	ADA_TIME0	<- seq(from = 10,to = 12,by = 1)	# IDs with ADA never present
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
	FINAL_ALB3 <- c(1,4,7,10)	# ID's with FINAL_ALB == 3
	FINAL_ALB4 <- c(2,5,8,11)	# ID's with FINAL_ALB == 4
	FINAL_ALB5 <- c(3,6,9,12)	# ID's with FINAL_ALB == 5
	# Fill in final albumin values based on ID
		cov.data$FINAL_ALB <- NA	#Add a FINAL_ALB column
		cov.data$FINAL_ALB[cov.data$ID %in% FINAL_ALB3] <- 3	# FINAL_ALB == 3
		cov.data$FINAL_ALB[cov.data$ID %in% FINAL_ALB4] <- 4	# FINAL_ALB == 4
		cov.data$FINAL_ALB[cov.data$ID %in% FINAL_ALB5] <- 5	# FINAL_ALB == 5

# ------------------------------------------------------------------------------
# Data frame of individual characteristics
	pop.data <- lapply(cov.data,rep.int,times = length(TIME))
	pop.data <- as.data.frame(pop.data)
	pop.data <- pop.data[with(pop.data, order(pop.data$SIM,pop.data$ID)), ]
	pop.data$TIME <- TIME

# Calculate albumin concentrations for each individual for all time-points
# A linear function containing the baseline albumin (BASE_ALB) and their last albumin (FINAL_ALB)
	pop.data <- ddply(pop.data, .(SIM,ID), albumin.function)

# Flag if ADA are present for each individual for all time-points
# This assumes that once a person develops ADA, they stay with ADA
	pop.data <- ddply(pop.data, .(SIM,ID), ada.function)

# Create a data frame of covariate values for each individual at key time-points
# i.e., baseline (0), day 98, 210, 378 and 546
	cov.time.data <- subset(pop.data,select = c(ID,SIM,TIME,ALB,ADA))
	cov.time.data <- cov.time.data[cov.time.data$TIME %in% c(0,98,210,378,546),]
