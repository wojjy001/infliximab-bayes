# Create test population to be input for Bayesian estimation of individual PK parameters
# Same population to be input into NONMEM and R/mrgsolve/optim
# The population will just be administered label dosing, i.e., induction phase: 5 mg/kg at day 0, 14 and 42, and maintenance phase: 5 mg/kg every 56 days
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))

# Set working directory
 	# work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/NONMEM_Check/"	# Mac directory
	work.dir <- "E:/Wojciechowski/infliximab-bayes/NONMEM_Check/"	# Server directory
	setwd(work.dir)

# Set seed for reproducible results
 	set.seed(123456)

# ------------------------------------------------------------------------------
# Set up test population
 	n <- 100	# Number of individuals to be simulated
	ID <- 1:n	# ID sequence

# Define event times and time sequence
	time <- seq(from = 0,to = 546,by = 7)	# General time sequence for simulation
	inf.times <- c(0,14,42,seq(from = 98,to = 546,by = 56))	# Infusion times (induction and maintenance phase dosing)

# Define population's covariate and random effect values
 	ALB <- rlnorm(n,meanlog = log(4),sdlog = 0.09)	# Albumin, g/dL
	plot(hist(ALB))
	WT <- rlnorm(n,meanlog = log(70),sdlog = 0.09)	# Weight, kg
	plot(hist(WT))

	ETA1 <- rnorm(n,mean = 0,sd = sqrt(0.106929))	# Random effect for CL
	ETA2 <- rnorm(n,mean = 0,sd = sqrt(0.0225))	# Random effect for V1
	ETA3 <- rnorm(n,mean = 0,sd = sqrt(1.21))	# Random effect for Q
	ETA4 <- rnorm(n,mean = 0,sd = sqrt(0.638401))	# Random effect for V2

#	Calculate the presence of ADA or not at each time depending on ADA_TIME
	ADA_TIME <- round(runif(n,min = 98,max = 600),digits = 0)	# Random onset in anti-drug antibodies (ADA)
	ada.data <- data.frame(ID,ETA1,ETA2,ETA3,ETA4,ALB,WT,ADA_TIME)

	ADA.onset.function <- function(ada.data) {
		ADA.rate <- c(0,1,1)
		ADA.times <- c(0,ada.data$ADA_TIME[1],600)
		ADA.onset <- approxfun(ADA.times,ADA.rate,method = "const")
		ADA <- ADA.onset(time)

		ada.data.rep <- lapply(ada.data,rep.int,times = length(time))
		ada.data.rep <- as.data.frame(ada.data.rep)	# Convert list to data frame
		ada.data.rep$ADA <- ADA # Add new ADA column
		ada.data.rep$time <- time	# Add time column
		ada.data.rep <- ada.data.rep[,c(1,10,2:7,9)]	# Remove the old "ADA_TIME" column
		ada.data.rep	# Return the data frame
	}
	input.data <- ddply(ada.data, .(ID), ADA.onset.function)

# Add pre-simulated residual error
 	input.data$ERRPRO <- rnorm(n*length(time),mean = 0,sd = sqrt(0.175561))

# ------------------------------------------------------------------------------
# Add dosing information to input data frame
	input.data$amt <- 0	# Initially make all values for amt = 0
	input.data$cmt <- 1	# Dosing into the central compartment
	input.data$evid <- 0	# Initially make only observations, no dosing
	input.data$rate <- 0	# Initially no rate specifications as dosing is 0

	input.data$amt[input.data$time %in% inf.times] <- input.data$WT[input.data$time %in% inf.times]*5	# 5 mg/kg dosing
	input.data$evid[input.data$time %in% inf.times] <- 1	# Specify dose is being administered
	input.data$rate[input.data$time %in% inf.times] <- -2	# Infusion duration is specified in the mrgsolve model code file

# Save data frame to .csv
	write.csv(input.data,file = paste0(work.dir,"mrgsolve_simulation_input.csv"),na = ".",quote = F,row.names = F)

# ------------------------------------------------------------------------------
# Make data frame NONMEM compatible for simulation
	nonmem.input.data <- input.data
	names(nonmem.input.data)[c(1,2,11:14)] <- c("CID","TIME","AMT","CMT","EVID","RATE")
	nonmem.input.data$DV <- "."
# NONMEM requires separate lines for observation and dosing event
	extra.obs.data <- nonmem.input.data[nonmem.input.data$TIME %in% inf.times,]
	extra.obs.data$AMT <- 0	# Value for no dose
	extra.obs.data$EVID <- 0	# Value for observation
	extra.obs.data$RATE <- 0
# Combine with input data
	nonmem.input.data <- rbind(nonmem.input.data,extra.obs.data)
	nonmem.input.data <- nonmem.input.data[with(nonmem.input.data, order(nonmem.input.data$CID,nonmem.input.data$TIME,-nonmem.input.data$EVID)),]	
	write.csv(nonmem.input.data,file = paste0(work.dir,"nonmem_simulation_input.csv"),na = ".",quote = F,row.names = F)
