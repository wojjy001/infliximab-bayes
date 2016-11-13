# in silico infliximab dosing project
# Script for containing universal functions
# ------------------------------------------------------------------------------
# Create and set working directory
# Specific for the simulation
	n <- 9	# Number of seed individuals (where each seed individual has a different set of covariate values)
	nsim <- 10	# Number of simulations of the seed individuals to perform
	# Set seed for reproducible results
		seed <- round(runif(1,min = 000001,max = 999999),digits = 0)
		# seed <- 647914
		set.seed(seed)
	# sim.name <- paste("SUCCESS_SIM",nsim,"_IND",n,"_seed",seed,sep = "")	# Simulation folder's name
	sim.name <- paste("SIM",nsim,"_IND",n,"_seed",seed,sep = "")	# Simulation folder's name
	# sim.output.dir <- paste0("D:/Moved-Infliximab-Output/",sim.name,"/")	# Simulation directory for Windows
	sim.output.dir <- paste0("E:/Wojciechowski/Moved-Infliximab-Output/",sim.name,"/")	# Simulation directory for Server
	# sim.output.dir <- paste0("/Volumes/Prosecutor/PhD/InfliximabBayes/Moved-Infliximab-Output/",sim.name,"/")	# Simulation directory for Mac
	dir.create(file.path(sim.output.dir),showWarnings = FALSE) # Create simulation directory
	setwd(file.path(sim.output.dir))	#Set the working directory

# ------------------------------------------------------------------------------
# Load package libaries
	library(ggplot2)	# Plotting package
	library(grid)	# Plotting package
# Custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 14))

# ------------------------------------------------------------------------------
# Pre-defined universal objects
# Target trough concentration definitions
	trough.target <- 3	# Set the target trough concentration for dose optimisation
	trough.upper <- 5	# Set upper bound for trough concentrations

# Values for PPV (Population Parameter Variability), as SDs
	PPVCL <- 0.1345
	PPVV1 <- 0.4207
	PPVQ <- 0.8515
	PPVV2 <- 0.3240
# Value for RUV (Residual Unexplained Variability), as SD
	ERRPRO <- 0.3

# Define time sequences
	# Infusion times (0, 2, 6 weeks and then every 8 weeks) in days
		TIME1i <- c(0,14,42)
		TIMEi <- c(TIME1i,98,154,210,266,322,378,434,490,546)
	# Infusion duration (2 hours) in days
		INFD <- 2/24
	# Overall time sequence
		time.int <- 1	# Difference in simulation times
		TIME <- seq(from = 0,to = 602,by = time.int)
	# Object specifying beyond the TIME sequence
		END <- max(TIME)+100
	# Define the last time-point to be simulated
		last.time <- 602	# days
	# After the initiation phase, the first sample will be collected at day 98
		sample.times <- c(0,98)	# days
	# Initial dosing interval for the maintenance phase
		dose.int <- 56	# days
		next.dose.int <- 56	# days

# Set the dose for simulating the first intervals
	amt.init1 <- 5	# initial dose mg/kg
	amt.init2 <- 10	# initial dose mg/kg
# Set the min and max mg/kg doses for bayesian dosing
	amt.min <- 3	# minimum dose mg/kg
	amt.max <- 10	# maximum dose mg/kg

# ------------------------------------------------------------------------------
# Pre-defined universal functions
# Function for calculating changes in random effects
# A linear function containing the baseline ETA (BASE_ETA) and their last ETA (FINAL_ETA)
	eta.function <- function(input.data) {
		if (time.dep == 0) {
			input.data$ETA2 <- input.data$BASE_ETA2
			input.data$ETA3 <- input.data$BASE_ETA3
			input.data$ETA4 <- input.data$BASE_ETA4
		}
		if (time.dep == 1) {
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
		}
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
