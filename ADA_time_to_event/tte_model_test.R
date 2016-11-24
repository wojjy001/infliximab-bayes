# in silico infliximab dosing project
# Onset of ADA time to event model
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list=ls(all=TRUE))

# Load package libraries
	library(ggplot2)	# Plotting
	library(grid)	# Plotting
	library(plyr) # Split and rearrange data, ddply function
	library(dplyr)	# Split and rearrange data
	library(mrgsolve)	# Metrum Research Group differential equation solver

# Set working directory
	setwd("/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/ADA_time_to_event/")

# ------------------------------------------------------------------------------
# Define the model parameters and equations
# Using mrgsolve - differential equations
# This compiled model is used for simulating n individuals and their concentration-time profiles
	code <- '
	$INIT			// Initial conditions for compartments
						CENT = 0,	// Central for PK
						PERI = 0, // Peripheral for PK
						CMHAZ = 0,	// Cumulative hazard
						TAT = 0,	//	Time above target trough concentration
						AUT = 0,	// Area under the target trough concentration
						TBT = 0	// Time below target trough concentration

	$PARAM		// Population parameters
						POPCL = 0.294,
						POPV1 = 3.33,
						POPQ = 0.0719,
						POPV2 = 1.14,
						POPB0 = -2.11,
						POPB1 = - 7.07,
						POPET50 = 4.62,

						// Covariate effects
						WT_CL = 0.614,	// Effect of weight on clearance
						WT_V1 = 0.691,	// Effect of weight on V1
						WT_Q = 1.1,	// Effect of weight on Q
						WT_V2 = 0.59,	// Effect of weight on V2
						ALB_CL = -1.17,	// Effect of albumin on clearance
						ADA_CL = 0.257,	// Effect of anti-drug antibodies on clearance
						ALB_HAZ = -0.0479,	// Effect of albumin on hazard
						WT_HAZ = -0.0162,	// Effect of weight on hazard
						SL_DE = -2.93,	// Drug effect slope on hazard
						CP_HAZ = -2.65,	// Time-varying hazard

						// Covariate values for simulation
						WT = 70,	// Weight (kg)
						ALB = 4,	// Albumin (U/L)
						target = 3,	// Target trough concentration (mg/L)
						SIM = 0,	// Simulation identifier
						RANDOM = 0,	// Random number pulled from uniform distribution (0,1)

						// Presimulated PPV values
						ETA1 = 0,
						ETA2 = 0,
						ETA3 = 0,
						ETA4 = 0,
						ERRPRO = 0

	$OMEGA		// Inter-individual variability
						name = "BSV"
						block = FALSE
						labels = s(PPVCL,PPVV1,PPVQ,PPVV2)
						0.106929
						0.0225
						1.21
						0.638401

	$SIGMA		// Intra-individual variability
						block = FALSE
						labels = s(ERR_PRO)
						0.175561

	$MAIN			// Infusion duration
						D_CENT = 0.08333333;	// 2 hours

						// Set the relative time since last dose (RTLD)
						if (NEWIND <= 1) {
							double BASEALB = ALB;	// Initialise baseline albumin for each new individual
							double BASEWT = WT;	// Initialise baseline weight for each new individual
							double OLDCMHAZ = 0;	// Old cumulative hazard
						}

						// Covariate effects
						double NWT = WT/70;	// Normalised weight (kg)
						double NALB = ALB/4;	// Normalised albumin (U/L)
						double ADACOV = 1;	// No anti-drug antibodies
						if (ADA == 1) ADACOV = 1+ADA_CL;		// Anti-drug antibodies

						// Individual pharmacokinetic parameter values
						double CL = POPCL*pow(NWT,WT_CL)*pow(NALB,ALB_CL)*ADACOV*exp(ETA1);
						double V1 = POPV1*pow(NWT,WT_V1)*exp(ETA2);
						double Q = POPQ*pow(NWT,WT_Q)*exp(ETA3);
						double V2 = POPV2*pow(NWT,WT_V2)*exp(ETA4);

						// Micro-rate constants for infliximab pharmacokinetics
						double K10 = CL/V1;
						double K12 = Q/V1;
						double K21 = Q/V2;

						// Half-life
						double Thalf = log(2)/(0.5*((K10+K12+K21)-sqrt(pow(K10+K12+K21,2)-4*K10*K21)));

						// Parameters for time to formation of ADA
						double B0 = POPB0;	// Baseline hazard
						double B1 = POPB1;	// Hazard slope - linear with time
						double ET50 = exp(POPET50);	// Time to achieve 50% of hazard?
						double B2 = ALB_HAZ;
						double B3 = WT_HAZ;
						double SL = SL_DE;	// Drug effect slope on hazard
						double CP = 1/(1+exp(CP_HAZ));	// Plasma concentration threshold

	$ODE			// Differential equations for infliximab pharmacokinetics
						dxdt_CENT = -K12*CENT +K21*PERI -K10*CENT;
						dxdt_PERI = K12*CENT -K21*PERI;

						// Plasma concentration
						double C1 = CENT/V1;	// Plasma concentration of the central compartment

						// Time above, area under and time below target concentration
						dxdt_TAT = 0;
						dxdt_AUT = 0;
						dxdt_TBT = 0;
						if (SOLVERTIME > 0.08333333 & C1 > target) dxdt_TAT = 1;
						if (SOLVERTIME > 0.08333333 & C1 < target) {
							dxdt_AUT = target - C1;
							dxdt_TBT = 1;
						}

						// Equations for time to formation of ADA
						double PROP = 1;
						if (SOLVERTIME > 0.08333333) PROP = TAT/SOLVERTIME;	// Proportion of time above target
						double DE = 0; // Drug effect
						if (PROP > CP) DE = 1;
						double LH = B0+B1*SOLVERTIME/(ET50+SOLVERTIME) +B2*(BASEALB-4) +B3*(BASEWT-70) +SL*DE;	// log-hazard
						double HAZD = exp(LH);	// Hazard

						// Calculate survivor function
						dxdt_CMHAZ = HAZD;	// Cumulative hazard from time = 0
						double LASTCMHAZ = CMHAZ - OLDCMHAZ; // Cumulative hazard since last time
						OLDCMHAZ = CMHAZ;	// Resetting the old cumulative hazard
						double ST = exp(-CMHAZ);	// Survivor function

						// Determine presence of anti-drug antibodies
						double ADA = 0;
						if (RANDOM >= ST) ADA = 1;

	$TABLE		// Predicted concentrations
						table(IPRE) = CENT/V1;
						table(DV) = table(IPRE)*(1+ERRPRO);

	$CAPTURE	// Variables to output
						WT ADA ALB RANDOM CL V1 Q V2 ETA1 ETA2 ETA3 ETA4 PROP LASTCMHAZ OLDCMHAZ ST
	'
# Compile the model code
	mod <- mcode("popADATTE",code)
		# There is opportunity to simply update model parameters after the model code has been compiled

# ------------------------------------------------------------------------------
#	Create input data frame
	inf.times <- c(0,14,42,98,154,210,266,322,378,434,490)
	times <- seq(from = 0,to = 546,by = 7)
	input.data <- data.frame(
		ID = 1,
		time = times,
		amt = 175,
		evid = 1,
		cmt = 1,
		rate = -2,
		ALB = 4,
		RANDOM = runif(length(times),min = 0,max = 1)
	)
	# If not an infusion time, make amt, evid, and rate = 0
	 	input.data$amt[!c(input.data$time %in% inf.times)] <- 0
	 	input.data$evid[!c(input.data$time %in% inf.times)] <- 0
	 	input.data$rate[!c(input.data$time %in% inf.times)] <- 0

# Simulate
	sim.data <- mod %>% mrgsim(data = input.data) %>% as.tbl
	sim.data <- as.data.frame(sim.data)
	head(sim.data)

# ------------------------------------------------------------------------------
#	Plot results
# Infliximab concentrations
	plotobj1 <- NULL
	plotobj1 <- ggplot(sim.data)
	plotobj1 <- plotobj1 + geom_line(aes(x = time,y = IPRE),colour = "red")
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 3),linetype = "dashed")
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 5),linetype = "dashed")
	plotobj1 <- plotobj1 + scale_y_log10("Infliximab Concentration (mg/L)\n",breaks = c(0.001,0.01,0.1,1,10,100,1000))
	plotobj1 <- plotobj1 + scale_x_continuous("\nTime (days)")
	print(plotobj1)

#	Hazard and survival
 	plotobj2 <- NULL
	plotobj2 <- ggplot(sim.data)
	plotobj2 <- plotobj2 + geom_line(aes(x = time,y = ST),colour = "darkgreen")
	plotobj2

table(sim.data$ADA)
