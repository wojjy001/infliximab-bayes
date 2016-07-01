# Time-weighted Bayes project
# Script for containing universal functions
#-------------------------------------------------------------------------------
# Load package libaries
	library(ggplot2)	# Plotting package
	library(grid)	# Plotting package
	library(plyr)	# Split and rearrange data, ddply functions
	library(mrgsolve)	# Metrum Research Group differential equation solver for pharmacometrics

# Custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 14))

# Set seed for reproducible results
	set.seed(123456)

#-------------------------------------------------------------------------------
# Pre-defined universal objects
# Target trough concentration definitions
	trough.target <- 3	# Set the target trough concentration for dose optimisation
	trough.upper <- 5	# Set upper bound for trough concentrations

# Values for PPV (Population Parameter Variability), as SDs
	PPVCL <- 0.327
	PPVV1 <- 0.150
	PPVQ <- 1.10
	PPVV2 <- 0.799

# Define time sequences
	# Infusion times (0, 2, 6 weeks and then every 8 weeks) in days
		TIMEi <- c(0,14,42,98,154,210,266,322,378,434,490)
	# Infusion duration (2 hours) in days
		INFD <- 2/24
	# Time sequence for the different sampling intervals, days
		TIME1 <- seq(from = 0,to = 98,by = 14)
		# Times in TIME1 that are infusion times
			TIME1i <- TIME1[TIME1 %in% TIMEi]
		TIME2 <- seq(from = 98,to = 210,by = 14)
		# Times in TIME2 that are infusion times
			TIME2i <- TIME2[TIME2 %in% TIMEi]
		TIME3 <- seq(from = 210,to = 378,by = 14)
		# Times in TIME3 that are infusion times
			TIME3i <- TIME3[TIME3 %in% TIMEi]
		TIME4 <- seq(from = 378,to = 546,by = 14)
		# Times in TIME4 that are infusion times
			TIME4i <- TIME4[TIME4 %in% TIMEi]

	# Overall time sequence
		TIME <- unique(sort(c(TIME1,TIME2,TIME3,TIME4)))
	# Object specifying beyond the TIME sequence
		END <- max(TIME)+100

#-------------------------------------------------------------------------------
# Pre-defined universal functions
# Functions for calculating 95% prediction intervals
	CI95lo <- function(x) quantile(x,probs = 0.025)
	CI95hi <- function(x) quantile(x,probs = 0.975)

# Function for taking the last row of a given factor (commonly use for taking the last row of each individual)
	lastperID <- function(x) tail(x,1)

# ------------------------------------------------------------------------------
# Define the model parameters and equations
# Using mrgsolve - differential equations
# This compiled model is used for simulating n individuals and their concentration-time profiles
	code <- '
	$INIT			// Initial conditions for compartments
						CENT = 0,	// Central
						PERI = 0, // Peripheral
						AUT = 0	//Time below target compartment

	$PARAM		// Population parameters
						POPCL = 0.294,
						POPV1 = 3.33,
						POPQ = 0.0719,
						POPV2 = 1.14,

						// Covariate effects
						WT_CL = 0.614,	// Effect of weight on clearance
						WT_V1 = 0.691,	// Effect of weight on V1
						WT_Q = 1.1,	// Effect of weight on Q
						WT_V2 = 0.59,	// Effect of weight on V2
						ALB_CL = -1.17,	// Effect of albumin on clearance
						ADA_CL = 0.257,	// Effect of anti-drug antibodies on clearance

						// Covariate values for simulation
						WT = 70,	// Weight (kg)
						ALB = 4,	// Albumin
						ADA = 0,	// Anti-drug antibodies (0 = No, 1 = Yes)
						target = 3	// Target trough concentration (mg/L)
						SIM = 0	// Simulation identifier

						// Presimulated PPV values
						ETA1 = 0,
						ETA2 = 0,
						ETA3 = 0,
						ETA4 = 0

	$OMEGA		name = "BSV"
						block = FALSE
						labels = s(PPVCL,PPVV1,PPVQ,PPVV2)
						0.106929
						0.0225
						1.21
						0.638401

	$SIGMA		block = FALSE
						labels = s(ERRPRO)
						0.175561

	$MAIN			// Infusion duration
						D_CENT = 0.08333333;  // 2 hours

						// Covariate effects
						double ADACOV = 1;	// No anti-drug antibodies
						if (ADA == 1) ADACOV = 1+ADA_CL;		// Anti-drug antibodies

						// Individual parameter values
						double CL = POPCL*pow(WT/70,WT_CL)*pow(ALB/4,ALB_CL)*ADACOV*exp(ETA1);
						double V1 = POPV1*pow(WT/70,WT_V1)*exp(ETA2);
						double Q = POPQ*pow(WT/70,WT_Q)*exp(ETA3);
						double V2 = POPV2*pow(WT/70,WT_V2)*exp(ETA4);

	$ODE			// Differential equations
						dxdt_CENT = -Q/V1*CENT +Q/V2*PERI -CL/V1*CENT;
						dxdt_PERI = Q/V1*CENT -Q/V2*PERI;

						// Time below target
						double CP = CENT/V1;	// Plasma concentration of the central compartment
						dxdt_AUT = 0;
						if (CP < target) dxdt_AUT = 1;

	$TABLE		table(IPRE) = CENT/V1;
						table(DV) = table(IPRE)*exp(ERRPRO);

	$CAPTURE	SIM WT ADA ALB CL V1 Q V2 ETA1 ETA2 ETA3 ETA4
	'
# Compile the model code
	mod <- mcode("popINFLIX",code)
		# There is opportunity to simply update model parameters after the model code has been compiled
