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
	$INIT			// Initial conditions for PK compartments
						CENT = 0,
						PERI = 0

	$PARAM		// Population parameters
						POPCL = 0.294,
						POPV1 = 3.33,
						POPQ = 0.0719,
						POPV2 = 1.14,

						// Covariate effects
						WT_CL = 0.614,	// Effect of weight on clearance
						ALB_CL = -1.17,	// Effect of albumin on clearance
						ADA_CL = 0.257,	// Effect of anti-drug antibodies on clearance
						WT_V1 = 0.691,	// Effect of weight on V1
						WT_Q = 1.1,	// Effect of weight on Q
						WT_V2 = 0.59,	// Effect of weight on V2

						// Covariate values
						ALB = 4,	// Albumin, g/dL
						WT = 70,	// Weight, kg
						ADA = 0, // Anti-drug antibodies, 0 = no, 1 = yes

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
							// Anti-drug antibodies
							double ADACOV = 1;	// No anti-drug antibodies
							if (ADA == 1) ADACOV = 1+ADA_CL; // Anti-drug antibodies

						// Individual parameter values
						double CL = POPCL*pow(WT/70,WT_CL)*pow(ALB/4,ALB_CL)*ADACOV*exp(ETA1);
						double V1 = POPV1*pow(WT/70,WT_V1)*exp(ETA2);
						double Q = POPQ*pow(WT/70,WT_Q)*exp(ETA3);
						double V2 = POPV2*pow(WT/70,WT_V2)*exp(ETA4);

						// Micro-rate constants
						double K10 = CL/V1;
						double K12 = Q/V1;
						double K21 = Q/V2;

	$ODE			// Differential equations
						dxdt_CENT = -K12*CENT +K21*PERI -K10*CENT;
						dxdt_PERI = K12*CENT -K21*PERI;

	$TABLE		table(IPRE) = CENT/V1;
						table(DV) = table(IPRE)*(1+ERRPRO);

	$CAPTURE	WT ALB ADA CL V1 Q V2 ETA1 ETA2 ETA3 ETA4
	'
# Compile the model code
	mod <- mcode("popINFLIX",code)
		# There is opportunity to simply update model parameters after the model code has been compiled
