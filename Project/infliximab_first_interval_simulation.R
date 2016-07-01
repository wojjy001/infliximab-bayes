# Time-weighted Bayes project
# Script for simulating concentrations for the first intervals
# Everyone will receive 5 mg/kg doses for the first interval in all scenarios
#-------------------------------------------------------------------------------
# Source the other R scripts and execute
	# infliximab_population.R
	source(paste0(work.dir,"infliximab_population.R"))

#-------------------------------------------------------------------------------
# Set the dose for simulating the first intervals
	amt1 <- 5	# 5 mg/kg

# Time sequence for the first interval in mrgsolve format
	tgrid1 <- tgrid(head(TIME1,1),tail(TIME1,1),14)

# Create a data frame ready for mrgsolve simulation
	# Function for creating a data frame ready for mrgsolve simulation
		population.interval1.function <- function(input.data) {
			ID.number <- input.data$ID[1]	# Individual ID
			SIM.number <- input.data$SIM[1]	# Individual simulation number
			ALB <- input.data$ALB[input.data$TIME %in% TIME1i]	# Individual albumin
			ADA <- input.data$ADA[input.data$TIME %in% TIME1i]	# Individual ADA status
			ETA1 <- input.data$ETA1[input.data$TIME %in% TIME1i]
			ETA2 <- input.data$ETA2[input.data$TIME %in% TIME1i]
			ETA3 <- input.data$ETA3[input.data$TIME %in% TIME1i]
			ETA4 <- input.data$ETA4[input.data$TIME %in% TIME1i]
			input.conc.data1 <- data.frame(
				ID = ID.number,
				SIM = SIM.number,
				time = TIME1i,	# Infusion times
				ALB,	# Albumin
				ADA,	# Anti-drug antibodies
				ETA1,
				ETA2,
				ETA3,
				ETA4,
				amt = amt1*70,	# mg/kg dose
				evid = 1,	# Dosing event
				cmt = 1,	# Dose into the central compartment (compartment = 1)
				rate = -2	# Infusion duration is specific in the model file
			)
			# Make the amt given in the last time-point == 0
				input.conc.data1$amt[input.conc.data1$time == max(TIME1i)] <- 0
				input.conc.data1$evid[input.conc.data1$time == max(TIME1i)] <- 0
				input.conc.data1$rate[input.conc.data1$time == max(TIME1i)] <- 0
			# Return input.conc.data
				input.conc.data1
		}
	# Create population data frame ready for mrgsolve simulation
		input.conc.data1 <- ddply(pop.data, .(ID,SIM), population.interval1.function)

# Simulate concentration-time profiles for individuals in input.conc.data
	per.simulation.function1 <- function(input.data) {
		conc.data1 <- mod %>% data_set(input.data) %>% carry.out(SIM,amt) %>% mrgsim(tgrid = tgrid1)
		conc.data1 <- as.data.frame(conc.data1)
	}
	conc.data1 <- ddply(input.conc.data1, .(SIM), per.simulation.function1)
