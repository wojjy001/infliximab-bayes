# Time-weighted Bayes project
# Script for simulating concentrations for the first interval
# Everyone will receive 5 mg/kg doses for the first interval in all scenarios
# ------------------------------------------------------------------------------
# Source the other R scripts and execute
	work.dir <- "D:/infliximab-bayes/Project/"
	source(paste0(work.dir,"population.R"))
	source(paste0(work.dir,"model.R"))

# ------------------------------------------------------------------------------
# Set the dose for simulating the first intervals
	amt1 <- 5	# 5 mg/kg

# Create a data frame ready for mrgsolve simulation
	# Function for creating a data frame ready for mrgsolve simulation
		interval1.label <- function(input.data) {
			ID.number <- input.data$ID[1]	# Individual ID
			SIM.number <- input.data$SIM[1]	# Individual simulation number
			ALB <- input.data$ALB[input.data$TIME %in% TIME1]	# Individual albumin
			ADA <- input.data$ADA[input.data$TIME %in% TIME1]	# Individual ADA status
			ETA1 <- input.data$ETA1[input.data$TIME %in% TIME1]
			ETA2 <- input.data$ETA2[input.data$TIME %in% TIME1]
			ETA3 <- input.data$ETA3[input.data$TIME %in% TIME1]
			ETA4 <- input.data$ETA4[input.data$TIME %in% TIME1]
			ERRPRO <- input.data$ERRPRO[input.data$TIME %in% TIME1]
			input.conc.data <- data.frame(
				ID = ID.number,
				SIM = SIM.number,
				time = TIME1,	# Time points for simulation
				ALB,	# Albumin
				ADA,	# Anti-drug antibodies
				ETA1,
				ETA2,
				ETA3,
				ETA4,
				ERRPRO,
				amt = amt1*70,	# mg/kg dose
				evid = 1,	# Dosing event
				cmt = 1,	# Dose into the central compartment (compartment = 1)
				rate = -2	# Infusion duration is specific in the model file
			)
			# Make the amt given in the last time-point == 0
				input.conc.data$amt[!c(input.conc.data$time %in% TIME1i) | input.conc.data$time == max(TIME1i)] <- 0
				input.conc.data$evid[!c(input.conc.data$time %in% TIME1i) | input.conc.data$time == max(TIME1i)] <- 0
				input.conc.data$rate[!c(input.conc.data$time %in% TIME1i) | input.conc.data$time == max(TIME1i)] <- 0
			# Return input.conc.data
				input.conc.data
		}
	# Create population data frame ready for mrgsolve simulation
		input.conc.data <- ddply(pop.data, .(ID,SIM), interval1.label)

# Simulate concentration-time profiles for individuals in input.conc.data
	conc.data <- ddply(input.conc.data, .(SIM), conc.per.simulation)
	conc.data.x <- conc.data[conc.data$time < 98,]
