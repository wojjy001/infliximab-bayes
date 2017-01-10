# in silico infliximab dosing project
# Script for simulating concentrations for the first interval (i.e., initiation phase)
# Everyone will receive 5 mg/kg doses for the first interval in all scenarios
# ------------------------------------------------------------------------------
# Create a data frame ready for mrgsolve simulation
	# Function for creating a data frame ready for mrgsolve simulation
		first.int.function1 <- function(pop.data) {
			# Create input data frame ready for mrgsolve for simulation the initiation phase
				input.first.int.data <- data.frame(
					ID = pop.data$ID[1],	# Individual ID
					time = TIME,	# Time points for simulation
					BASE_WT = pop.data$BASE_WT[1],	# Baseline weight
					BASE_ALB = pop.data$BASE_ALB[1],	# Baseline albumin
					TIME_WT = pop.data$BASE_WT[1],	# Default weight at TIME
					TIME_ADA = 0,	# Default ADA status at TIME
					TIME_ALB = pop.data$BASE_ALB[1],	# Default albumin at TIME
					ADAr = pop.data$ADAr[1],
					ETA1 = pop.data$ETA1,	# Random effect for CL
					ETA2 = pop.data$ETA2,	# Random effect for V1
					ETA3 = pop.data$ETA3,	# Random effect for Q
					ETA4 = pop.data$ETA4,	# Random effect for V2
					ERRPRO = pop.data$ERRPRO,	# Residual error
					amt = amt.init1*BASE_WT,	# mg/kg dose
					evid = 1,	# Dosing event
					cmt = 1,	# Dose into the central compartment (compartment = 1)
					rate = -2	# Infusion duration is specified in the model file
				)
			# Make the amt in times that aren't infusion times == 0
				input.first.int.data$amt[!c(input.first.int.data$time %in% TIME1i)] <- 0
			# Specify the correct evid and rate for non-infusion times
				input.first.int.data$evid[!c(input.first.int.data$time %in% TIME1i)] <- 0
				input.first.int.data$rate[!c(input.first.int.data$time %in% TIME1i)] <- 0
			# Flag if we want covariates to change depending on concentrations during simulation
				input.first.int.data$FLAG <- time.dep
			# Simulate concentration-time profiles for individuals in input.conc.data
				first.int.data <- mod %>% mrgsim(data = input.first.int.data,carry.out = c("amt","ERRPRO")) %>% as.tbl
		}
# Simulate the initiation phase from "pop.data" using "first.int.function1" for each ID (individual) in each SIM (simulation group)
	first.int.data1 <- ddply(pop.data, .(ID,SIM), first.int.function1)
