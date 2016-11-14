# in silico infliximab dosing project
# Script for simulating characteristics of the study population
# Does not simulation concentrations
# ------------------------------------------------------------------------------
# Source the other R scripts and execute
	# work.dir <- "D:/infliximab-bayes/Project/"	# Windows directory
	work.dir <- "E:/Wojciechowski/infliximab-bayes/Project/"	# Server directory
	# work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/"	# Mac directory

# ------------------------------------------------------------------------------
# Define population's characteristics
# Only going to pre-specify weight as 70 kg and randomly generate ADA_TIME, BASE_ALB and FINAL_ALB
	ID <- seq(from = 1,to = n,by = 1)	# Sequence of individual IDs
	SIM <- sort(c(rep(seq(from = 1,to = nsim,by = 1),times = n)))	# Sequence of simulation identifiers

# Create a data frame - one row per individual of covariate and random effects
	ID.data <- data.frame(SIM,ID)
	cov.data <- ID.data	# Transfer to a data frame called "cov.data"

# Assign ID values to specific groups of covariate values
# WEIGHT
	BASE_WT40 <- c(1,2,3)	# ID's with baseline weight = 40 kg
	BASE_WT70 <- c(4,5,6)	# ID's with baseline weight = 70 kg
	BASE_WT100 <- c(7,8,9)	# ID's with baseline weight = 100 kg
	cov.data$BASE_WT[cov.data$ID %in% BASE_WT40] <- 40
	cov.data$BASE_WT[cov.data$ID %in% BASE_WT70] <- 70
	cov.data$BASE_WT[cov.data$ID %in% BASE_WT100] <- 100

# ALBUMIN
	BASE_ALB25 <- c(1,4,7)	# ID's with baseline albumin = 2.5 U/L
	BASE_ALB30 <- c(2,5,8)	# ID's with baseline albumin = 3 U/L
	BASE_ALB35 <- c(3,6,9)	# ID's with baseline albumin = 3.5 U/L
	cov.data$BASE_ALB[cov.data$ID %in% BASE_ALB25] <- 2.5
	cov.data$BASE_ALB[cov.data$ID %in% BASE_ALB30] <- 3
	cov.data$BASE_ALB[cov.data$ID %in% BASE_ALB35] <- 3.5

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

# Random number for ADA onset
	cov.data$ADAr <- runif(n*nsim,min = 0,max = 1)

# Random number for albumin Emax
	cov.data$ALBr <- rlnorm(n*nsim,meanlog = log(1),sdlog = 0.3)

# Random number for weight Emax
	cov.data$WTr <- rlnorm(n*nsim,meanlog = log(4),sdlog = 0.3)

# ------------------------------------------------------------------------------
# Data frame of individual characteristics
	pop.data <- lapply(cov.data,rep.int,times = length(TIME))
	pop.data <- as.data.frame(pop.data)
	pop.data <- pop.data[with(pop.data, order(pop.data$SIM,pop.data$ID)), ]
	pop.data$TIME <- TIME

# Simulating the same residual error
	pop.data$ERRPRO <- rnorm(length(pop.data$TIME),mean = 0,sd = ERRPRO)	# Proportional residual error
