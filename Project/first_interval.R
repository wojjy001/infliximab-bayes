# Time-weighted Bayes project
# Script for simulating concentrations for the first interval
# Everyone will receive 5 mg/kg doses for the first interval in all scenarios
# ------------------------------------------------------------------------------
# Source the other R scripts and execute
	work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/"
	source(paste0(work.dir,"population.R"))
	source(paste0(work.dir,"model.R"))

# ------------------------------------------------------------------------------
# Create a data frame ready for mrgsolve simulation
	# Function for creating a data frame ready for mrgsolve simulation
		first.int.function <- function(input.data) {
			ID.number <- input.data$ID[1]	# Individual ID
			SIM.number <- input.data$SIM[1]	# Individual simulation number
			ALB <- input.data$ALB	# Individual albumin
			ADA <- input.data$ADA	# Individual ADA status
			ETA1 <- input.data$ETA1
			ETA2 <- input.data$ETA2
			ETA3 <- input.data$ETA3
			ETA4 <- input.data$ETA4
			ERRPRO <- input.data$ERRPRO
			input.first.int.data <- data.frame(
				ID = ID.number,
				SIM = SIM.number,
				time = TIME,	# Time points for simulation
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
				input.first.int.data$amt[!c(input.first.int.data$time %in% TIME1i)] <- 0
				input.first.int.data$evid[!c(input.first.int.data$time %in% TIME1i)] <- 0
				input.first.int.data$rate[!c(input.first.int.data$time %in% TIME1i)] <- 0
			# Simulate concentration-time profiles for individuals in input.conc.data
				first.int.data <- mod %>% mrgsim(data = input.first.int.data,carry.out = c("amt","ERRPRO")) %>% as.tbl
		}
# Create population data frame ready for mrgsolve simulation
	first.int.data <- ddply(pop.data, .(ID,SIM), first.int.function, .parallel = TRUE)
