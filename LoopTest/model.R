# Time-weighted Bayes project
# Infliximab population model code
#-------------------------------------------------------------------------------
# Define the model parameters and equations
# Using mrgsolve - differential equations
# This compiled model is used for simulating n individuals and their concentration-time profiles
	code <- '
	$INIT			// Initial conditions for compartments
						CENT = 0,	// Central
						PERI = 0, // Peripheral
						AUT = 0,	// Area under target trough compartment
						TBT = 0	// Time below target trough compartment

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
						target = 3,	// Target trough concentration (mg/L)
						SIM = 0,	// Simulation identifier

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

						// Plasma concentration
						double CP = CENT/V1;	// Plasma concentration of the central compartment

						// Area below target and time below target
						dxdt_AUT = 0;
						dxdt_TBT = 0;
						if (SOLVERTIME > 0.08333333 & CP < target) {
							dxdt_AUT = target - CP;
							dxdt_TBT = 1;
						}

	$TABLE		table(IPRE) = CENT/V1;
						table(DV) = table(IPRE)*(1+ERRPRO);

	$CAPTURE	WT ADA ALB CL V1 Q V2 ETA1 ETA2 ETA3 ETA4
	'
# Compile the model code
	mod <- mcode("popINFLIX",code)
		# There is opportunity to simply update model parameters after the model code has been compiled
