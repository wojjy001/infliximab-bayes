# in silico infliximab dosing project
# Infliximab population model code
# ------------------------------------------------------------------------------
# Define the model parameters and equations
# Using mrgsolve - differential equations
# This compiled model is used for simulating n individuals and their concentration-time profiles
	code <- '
	$SET			atol = 1e-8, rtol = 1e-8
						maxsteps = 100000

	$INIT			// Initial conditions for PK compartments
						CENT = 0,
						PERI = 0,
						AUT = 0,
						TUT = 0

	$CMT			// Specify covariate compartments
						ALB,	// Albumin, U/L
						WT	// Weight, kg

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
						BASE_WT = 70,	// Waseline weight (kg)
						BASE_ALB = 3,	// Baseline albumin at treatment initiation
						TIME_WT = 70,	// Time-dependent weight (kg)
						TIME_ADA = 0,	// Time-dependent ADA status
						TIME_ALB = 3, // Time-dependent albumin (U/L)
						target = 3,	// Target trough concentration (mg/L)
						SIM = 0,	// Simulation identifier
						FLAG = 0,	// Scenario identifier

						// Presimulated PPV values
						ETA1 = 0,
						ETA2 = 0,
						ETA3 = 0,
						ETA4 = 0,
						ERRPRO = 0

	$OMEGA		name = "BSV"
						block = FALSE
						labels = s(PPVCL,PPVV1,PPVQ,PPVV2)
						0.106929
						0.0225
						1.21
						0.638401

	$SIGMA		block = FALSE
						labels = s(ERR_PRO)
						0.175561

	$MAIN			// Infusion duration
						D_CENT = 0.08333333;  // 2 hours

						// Compartment initial conditions for covariates
						ALB_0 = BASE_ALB;
						WT_0 = BASE_WT;

						// Covariate effects
							// Anti-drug antibodies
							double ADACOV = 1;	// No anti-drug antibodies
							if (ADA == 1) ADACOV = 1+ADA_CL; // Anti-drug antibodies

						if (FLAG == 1) {	// For Bayesian scenarios
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

						// Half-life
						double Thalf = log(2)/(0.5*((K10+K12+K21)-sqrt(pow(K10+K12+K21,2)-4*K10*K21)));

						// Previous time under target concentration
						double prevTUT = TUT;

	$ODE			// Differential equations
						dxdt_CENT = -K12*CENT +K21*PERI -K10*CENT;
						dxdt_PERI = K12*CENT -K21*PERI;

						// Plasma concentration
						double CP = CENT/V1;	// Plasma concentration of the central compartment

						// Area below target and time under/above target
						dxdt_AUT = 0;
						dxdt_TUT = 0;
						if (SOLVERTIME > 0.08333333 & CP < target) {
							dxdt_AUT = target - CP;
							dxdt_TUT = 1;
						}

						// Proportion of time under target
						double pTUT = 0;
						if (SOLVERTIME > 0.08333333) pTUT = TUT/SOLVERTIME;
						double TUTdiff = TUT - prevTUT;
						double pTUTdiff = TUTdiff/7;

						// Albumin
						// dxdt_ALB = 0;
						// dxdt_WT = 0;
						// if (SOLVERTIME > 98) {
							dxdt_ALB = ALB*0.001;
							if (pTUTdiff > 0.05 & pTUTdiff <= 0.1) dxdt_ALB = 0;
							if (pTUTdiff > 0.1) dxdt_ALB = ALB*-0.001;
							// Weight
							dxdt_WT = WT*0.0005;
							if (pTUTdiff > 0.05 & pTUTdiff <= 0.1) dxdt_WT = 0;
							if (pTUTdiff > 0.1) dxdt_WT = WT*-0.0005;
						//}

						if (FLAG == 0) {	// Simulation when new dose administered
							// Limits on albumin
							double ALBCOV = ALB;
							if (ALBCOV < 2) ALBCOV = 2;
							if (ALBCOV > 6) ALBCOV = 6;
							// Limits on weight
							double WTCOV = WT;
							if (WTCOV < 35) WTCOV = 35;
							if (WTCOV > 110) WTCOV = 110;
						}

						// Anti-drug antibodies
						double ADA = 0;
						if (pTUT > 0.1 & SOLVERTIME > 0.08333333) ADA = 1;

	$TABLE		table(IPRE) = CENT/V1;
						table(DV) = table(IPRE)*(1+ERRPRO);

	$CAPTURE	BASE_WT BASE_ALB ADA CL ALBCOV WTCOV V1 Q V2 ETA1 ETA2 ETA3 ETA4 Thalf pTUT pTUTdiff
	'
# Compile the model code
	mod <- mcode("popINFLIX",code)
		# There is opportunity to simply update model parameters after the model code has been compiled
