# in silico infliximab dosing project
# Script for containing universal functions
# ------------------------------------------------------------------------------
# Create and set working directory
# Specific for the simulation
	n <- 20	# Number of seed individuals (where each seed individual has a different set of covariate values)
	nsim <- 100	# Number of simulations of the seed individuals to perform
	sim.name <- paste("SIM",nsim,"_IND",n,sep = "")	# Simulation folder's name
	sim.output.dir <- paste0("D:/Moved-Infliximab-Output/",sim.name,"/")	# Simulation directory
	dir.create(file.path(sim.output.dir),showWarnings = FALSE) # Create simulation directory
	setwd(file.path(sim.output.dir))	#Set the working directory

# ------------------------------------------------------------------------------
# Load package libaries
	library(ggplot2)	# Plotting package
	library(grid)	# Plotting package
	library(plyr)	# Split and rearrange data, ddply functions
	library(dplyr)
	library(mrgsolve)	# Metrum Research Group differential equation solver for pharmacometrics
	library(compiler)	# Compile repeatedly called functions
	library(numDeriv)
# Custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 14))
# Set seed for reproducible results
	set.seed(123456)

# ------------------------------------------------------------------------------
# Pre-defined universal objects
# Target trough concentration definitions
	trough.target <- 3	# Set the target trough concentration for dose optimisation
	trough.upper <- 5	# Set upper bound for trough concentrations

# Values for PPV (Population Parameter Variability), as SDs
	PPVCL <- 0.327
	PPVV1 <- 0.150
	PPVQ <- 1.10
	PPVV2 <- 0.799
# Value for RUV (Residual Unexplained Variability), as SD
	ERRPRO <- 0.419

# Define time sequences
	# Infusion times (0, 2, 6 weeks and then every 8 weeks) in days
		TIME1i <- c(0,14,42)
		TIMEi <- c(TIME1i,98,154,210,266,322,378,434,490,546)
	# Infusion duration (2 hours) in days
		INFD <- 2/24
	# Overall time sequence
		TIME <- seq(from = 0,to = 600,by = 1)
	# Object specifying beyond the TIME sequence
		END <- max(TIME)+100
	# Define the last time-point to be simulated
		last.time <- 546	# days
	# After the initiation phase, the first sample will be collected at day 98
		sample.times <- c(0,98)	# days
	# Initial dosing interval for the maintenance phase
		dose.int <- 56	# days
		next.dose.int <- 56	# days

# Set the dose for simulating the first intervals
	amt1 <- 5	# 5 mg/kg

# ------------------------------------------------------------------------------
# Pre-defined universal functions
# Function for calculating albumin concentrations for each individual for all time-points
# A linear function containing the baseline albumin (BASE_ALB) and their last albumin (FINAL_ALB)
	albumin.function <- function(input.data) {
		TIMEalb <- c(min(input.data$TIME),max(input.data$TIME))
		RATEalb <- c(head(input.data$BASE_ALB,1),head(input.data$FINAL_ALB,1))
		step.alb <- approxfun(TIMEalb,RATEalb,method = "linear")	# Linear function
		input.data$ALB <- step.alb(input.data$TIME)*(1+AMP_ALB1*sin(2*pi*FREQ_ALB1*input.data$TIME+PHASE_ALB1)+AMP_ALB2*sin(2*pi*FREQ_ALB2*input.data$TIME+PHASE_ALB2)+AMP_ALB3*sin(2*pi*FREQ_ALB3*input.data$TIME+PHASE_ALB3))	# Apply function to every time-point
		as.data.frame(input.data)
	}

# Function for flagging if ADA are present for each individual for all time-points
# This assumes that once a person develops ADA, they stay with ADA
	ada.function <- function(input.data) {
		TIMEada <- c(min(input.data$TIME),input.data$ADA_TIME[1],END)	# Specify times when ADA changes
		RATEada <- c(0,1,1)	# Specify the values for it to change to
		step.ada <- approxfun(TIMEada,RATEada,method = "const")	# Step function
		input.data$ADA <- step.ada(input.data$TIME)	# Apply function to every time-point
		as.data.frame(input.data)
	}

# Function for calculating changes in random effects
# A linear function containing the baseline ETA (BASE_ETA) and their last ETA (FINAL_ETA)
	eta.function <- function(input.data) {
		# ETA2
			TIMEeta2 <- c(min(input.data$TIME),max(input.data$TIME))
			RATEeta2 <- c(head(input.data$BASE_ETA2,1),head(input.data$FINAL_ETA2,1))
			step.eta2 <- approxfun(TIMEeta2,RATEeta2,method = "linear")	# Linear function
			input.data$ETA2 <- step.eta2(input.data$TIME)	# Apply function to every time-point
		# ETA3
			TIMEeta3 <- c(min(input.data$TIME),max(input.data$TIME))
			RATEeta3 <- c(head(input.data$BASE_ETA3,1),head(input.data$FINAL_ETA3,1))
			step.eta3 <- approxfun(TIMEeta3,RATEeta3,method = "linear")	# Linear function
			input.data$ETA3 <- step.eta3(input.data$TIME)	# Apply function to every time-point
		# ETA4
			TIMEeta4 <- c(min(input.data$TIME),max(input.data$TIME))
			RATEeta4 <- c(head(input.data$BASE_ETA4,1),head(input.data$FINAL_ETA4,1))
			step.eta4 <- approxfun(TIMEeta4,RATEeta4,method = "linear")	# Linear function
			input.data$ETA4 <- step.eta4(input.data$TIME)	# Apply function to every time-point
		# If SIM = 0, i.e, population typical patient, make ETAs equal zero
			if (input.data$SIM[1] == 0) {
				input.data$ETA1 <- 0
				input.data$ETA2 <- 0
				input.data$ETA3 <- 0
				input.data$ETA4 <- 0
			}
		# Return the resulting data frame
		input.data
	}
