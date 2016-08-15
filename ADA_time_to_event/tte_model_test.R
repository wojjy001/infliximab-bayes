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
						CENT1 = 0,	// Central for PK
						PERI1 = 0, // Peripheral for PK
						CMHAZ = 0,	// Cumulative hazard
						HZLAST = 0,	//
						CENT2 = 0,	// Central for PK
						PERI2 = 0,	// Peripheral for PK
						TA1 = 0,	//	Time above target trough concentration
						TA2 = 0, 	//	Time above target trough concentration
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
						ADAAMT_CL = 0.658,	// Effect of anti-drug antibodies on clearance
						LINES_CL = 0.0638	// Effect of cell line on clearance
						ALB_HAZ = -0.0479,	// Effect of albumin on hazard
						WT_HAZ = -0.0162,	// Effect of weight on hazard
						LE0_HAZ = 0, // Effect of cell line 0 on hazard
						LE2_HAZ = 0, // Effect of cell line 2 on hazard
						SL_DE = -2.93,	// Drug effect slope on hazard
						CP_HAZ = -2.65,	// Effect of plasma concentration (CP) on hazard

						// Covariate values for simulation
						AGE = 40,	// Age (years)
						WT = 70,	// Weight (kg)
						ALB = 4,	// Albumin (U/L)
						CRP = 10,	// C-reactive protein (g/L)
						LINE = 1,	// Cell line (0, 1 or 2 where 1 is standard)
						ADAAMT = 100,	// Anti-drug antibodies titre
						DISDUR = 15,	// Disease duration (days)
						target = 3,	// Target trough concentration (mg/L)
						SIM = 0,	// Simulation identifier

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
						D_CENT1 = 0.08333333;  // 2 hours

						// Set the relative time since last dose (RTLD)
						if (NEWIND <= 1) {
							double BASEALB = ALB;	// Initialise baseline albumin for each new individual
							double BASEWT = WT;	// Initialise baseline weight for each new individual
						}

						// Covariate effects
						double NWT = WT/70;	// Normalised weight (kg)
						double NAGE = AGE/40;	// Normalised age (years)
						double NALB = ALB/4;	// Normalised albumin (U/L)
						double NCRP = CRP/10;	// Normalised C-reactive protein (g/L)
						double NDISDUR = DISDUR/15;	// Normalised disease duration (days)
						double NADAAMT = ADAAMT/100;	// Normalised anti-drug antibody (ADA) titre
						double LINES = exp(LINES_CL*(LINE-1));	// Line standard is 1 (0, 1 or 2)
						double ADAAMTCOV = LINES*(1+pow(NADAAMT,ADAAMT_CL));	// Effect of ADA titre on CL

						// Individual pharmacokinetic parameter values
						double CL = POPCL*pow(NWT,WT_CL)*pow(NALB,ALB_CL)*ADAAMTCOV*exp(ETA1);
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
						double LE = 0; // Effect of cell line 1 on hazard
						if (LINE == 0) LE = LE0_HAZ;	// Effect of cell line 0 on hazard
						if (LINE == 2) LE = LE2_HAZ;	// Effect of cell line 2 on hazard
						double SL = SL_DE;	// Drug effect slope on hazard
						double CP = 1/(1+exp(CP_HAZ));

	$ODE			// Differential equations for infliximab pharmacokinetics
						dxdt_CENT1 = -K12*CENT1 +K21*PERI1 -K10*CENT1;
						dxdt_PERI1 = K12*CENT1 -K21*PERI1;
						dxdt_CENT2 = -K12*CENT2 +K21*PERI2 -K10*CENT2;
						dxdt_PERI2 = K12*CENT2 -K21*PERI2;

						// Plasma concentration
						double C1 = CENT1/V1;	// Plasma concentration of the central compartment
						double C5 = CENT2/V1;	// Plasma concentration of the central compartment

						// Equations for time to formation of ADA
						if (SOLVERTIME > 0) double PROP = TA1/SOLVERTIME;
						double DE = 0; // Drug effect
						if (PROP > CP) DE = 1;
						double LH = B0+B1*SOLVERTIME/(ET50+SOLVERTIME) +B2*(BASEALB-40) +B3*(BASEWT-70) +LE +SL*DE;	// log-hazard
						double HAZD = exp(LH);	// Hazard

						dxdt_CMHAZ = HAZD;
						dxdt_HZLAST = HAZD;
						dxdt_TA1 = 0;
						if (C1 > target) dxdt_TA1 = 1;
						dxdt_TA2 = 0;
						if (C5 > target) dxdt_TA2 = 1;

						// Area and time below target concentration
						dxdt_AUT = 0;
						dxdt_TBT = 0;
						if (SOLVERTIME > 0.08333333 & C1 < target) {
							dxdt_AUT = target - C1;
							dxdt_TBT = 1;
						}

	$TABLE		// Predicted concentrations
						table(IPRE) = CENT1/V1;
						table(DV) = table(IPRE)*(1+ERRPRO);

	$CAPTURE	// Variables to output
						WT ADAAMT ALB CL V1 Q V2 ETA1 ETA2 ETA3 ETA4 Thalf PROP CP
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
		amt = 350,
		evid = 1,
		cmt = 1,
		rate = -2
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

#	Hazard
 	plotobj2 <- NULL
	plotobj2 <- ggplot(sim.data)
	plotobj2 <- plotobj2 + geom_line(aes(x = time,y = CMHAZ),colour = "blue")
	plotobj2 <- plotobj2 + geom_line(aes(x = time,y = PROP),colour = "red")
	plotobj2
