# in silico infliximab dosing project
# Infliximab population model code
# ------------------------------------------------------------------------------
# Define the model parameters and equations
# Using mrgsolve - differential equations
	library(plyr)	# Split and rearrange data, ddply functions
	library(dplyr)
	library(mrgsolve)	# Metrum Research Group differential equation solver for pharmacometrics
# This compiled model is used for simulating n individuals and their concentration-time profiles
	code <- '
	$SET			atol = 1e-8, rtol = 1e-8
						maxsteps = 100000

	$INIT			// Initial conditions for PK compartments
						CENT = 0,
						PERI = 0,
						TUT = 0

	$CMT			// Specify covariate compartments
						ALB,	// Albumin, g/dL
						WT	// Weight, kg

	$PARAM		// Population parameters
						POPCL = 0.381,
						POPV1 = 2.37,
						POPQ = 0.122,
						POPV2 = 0.604,

						// Covariate effects
						WT_CL = 0.612,	// Effect of weight on clearance
						WT_V1 = 0.696,	// Effect of weight on V1
						WT_Q = 1.15,	// Effect of weight on Q
						WT_V2 = 0.604,	// Effect of weight on V2
						ALB_CL = -1.39,	// Effect of albumin on clearance
						ADA_CL = 1.59,	// Effect of anti-drug antibodies on clearance

						// Covariate values for simulation
						BASE_WT = 70,	// Baseline weight (kg)
						BASE_ALB = 3,	// Baseline albumin at treatment initiation
						TIME_WT = 70,	// Time-dependent weight (kg)
						TIME_ADA = 0,	// Time-dependent ADA status
						TIME_ALB = 3, // Time-dependent albumin (g/dL)
						ADAr = 0,	// ADA random number
						target = 3,	// Target trough concentration (mg/L)
						SIM = 0,	// Simulation identifier
						FLAG = 0,	// Time-dependence scenario identifier

						// Presimulated PPV values
						ETA1 = 0,
						ETA2 = 0,
						ETA3 = 0,
						ETA4 = 0,
						ERRPRO = 0

	$OMEGA		name = "BSV"
						block = FALSE
						labels = s(PPVCL,PPVV1,PPVQ,PPVV2)
						0.01809025
						0.1769885
						0.7250523
						0.104976

	$SIGMA		block = FALSE
						labels = s(ERR_PRO)
						0.09

	$MAIN			// Infusion duration
						D_CENT = 0.08333333;  // 2 hours

						// Compartment initial conditions for covariates
						ALB_0 = BASE_ALB;
						WT_0 = BASE_WT;

						// Covariate effects
							// Anti-drug antibodies
							double ADACOV = 1;	// No anti-drug antibodies
							if (ADA == 1) ADACOV = 1+ADA_CL; // Anti-drug antibodies

						if (FLAG == 0) {	// For Bayesian scenarios or non-time dependent scenarios
							ALBCOV = TIME_ALB;
							WTCOV = TIME_WT;
						}

						// Individual parameter values
						double CL = POPCL*pow(WTCOV/70,WT_CL)*pow(ALBCOV/4,ALB_CL)*ADACOV*exp(ETA1);
						double V1 = POPV1*pow(WTCOV/70,WT_V1)*exp(ETA2);
						double Q = POPQ*pow(WTCOV/70,WT_Q)*exp(ETA3);
						double V2 = POPV2*pow(WTCOV/70,WT_V2)*exp(ETA4);

						// Micro-rate constants
						double K10 = CL/V1;
						double K12 = Q/V1;
						double K21 = Q/V2;

						// Previous time under target concentration
						double prevTUT = TUT;

	$ODE			// Differential equations
						dxdt_CENT = -K12*CENT +K21*PERI -K10*CENT;
						dxdt_PERI = K12*CENT -K21*PERI;

						// Plasma concentration
						double CP = CENT/V1;	// Plasma concentration of the central compartment

						// Time under target
						dxdt_TUT = 0;
						if (SOLVERTIME > 0.08333333 & CP < target) dxdt_TUT = 1;

						// Proportion of time under target
						double pTUT = 0;
						if (SOLVERTIME > 0.08333333) {
							pTUT = TUT/SOLVERTIME;
						}
						double TUTdiff = TUT - prevTUT;
						double pTUTdiff = TUTdiff;

						// Albumin
						dxdt_ALB = 150*2.5/pow(150+SOLVERTIME,2);	// First derivative of Emax equation, Emax = 2.5 U/L, ET50 = 150 days
						if (pTUTdiff > 0.05 & pTUTdiff <= 0.1) dxdt_ALB = 0;
						if (pTUTdiff > 0.1) dxdt_ALB = 150*-2.5/pow(150+SOLVERTIME,2);  // First derivative of Emax equation, Emax = -2.5 U/L, ET50 = 150 days
						// Weight
						dxdt_WT = 150*8/pow(150+SOLVERTIME,2);	// First derivative of Emax equation, Emax = 8 kg, ET50 = 150 days
						if (pTUTdiff > 0.05 & pTUTdiff <= 0.1) dxdt_WT = 0;
						if (pTUTdiff > 0.1) dxdt_WT = 150*-8/pow(150+SOLVERTIME,2);	// First derivative of Emax equation, Emax = -8 kg, ET50 = 150 days

						if (FLAG == 1) {	// Simulation when new dose administered
							// Limits on albumin
							double ALBCOV = ALB;
							if (ALBCOV < 2) ALBCOV = 2;
							if (ALBCOV > 6) ALBCOV = 6;
							// Limits on weight
							double WTCOV = WT;
							if (WTCOV < 35) WTCOV = 35;
							if (WTCOV > 110) WTCOV = 110;
						}

						// Probability of developing anti-drug antibodies
						double ADA = 0;
						double ADAp = 0;
						if (FLAG == 1) {
							ADAp = pow(TUT,5)/(pow(19,5)+pow(TUT,5));
							if (ADAr < ADAp) ADA = 1;
						}

	$TABLE		table(IPRE) = CENT/V1;
						table(DV) = table(IPRE)*(1+ERRPRO);

	$CAPTURE	BASE_WT WTCOV BASE_ALB ALBCOV ADA ADAp ADAr CL V1 Q V2 ETA1 ETA2 ETA3 ETA4 pTUT pTUTdiff
	'
# Compile the model code
	mod <- mcode("popINFLIX",code)
		# There is opportunity to simply update model parameters after the model code has been compiled
