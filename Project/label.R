# Time-weighted Bayes project
# Script for simulating concentrations for the second, third and fourth intervals
# Everyone will receive 5 mg/kg dose
# ------------------------------------------------------------------------------
# Create a data frame ready for mrgsolve simulation
	# Function for creating a data frame ready for mrgsolve simulation
		all.interval.label <- function(input.data) {
			ID.number <- input.data$ID[1]	# Individual ID
			SIM.number <- input.data$SIM[1]	# Individual simulation number
			ALB <- input.data$ALB	# Individual albumin
			ADA <- input.data$ADA	# Individual ADA status
			ETA1 <- input.data$ETA1
			ETA2 <- input.data$ETA2
			ETA3 <- input.data$ETA3
			ETA4 <- input.data$ETA4
			ERRPRO <- input.data$ERRPRO
			input.label.data <- data.frame(
				ID = ID.number,
				SIM = SIM.number,
				time = TIME,	# Time points for simulation
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
				input.label.data$amt[!c(input.label.data$time %in% TIMEi)] <- 0
				input.label.data$evid[!c(input.label.data$time %in% TIMEi)] <- 0
				input.label.data$rate[!c(input.label.data$time %in% TIMEi)] <- 0
			# Return input.label.data
				input.label.data
		}
	# Create population data frame ready for mrgsolve simulation
		input.label.data <- ddply(pop.data, .(ID,SIM), all.interval.label)

# Simulate concentration-time profiles for individuals in input.label.data
	label.data <- ddply(input.label.data, .(SIM), conc.per.simulation)
	label.data <- label.data[with(label.data, order(label.data$ID,label.data$SIM)), ]	# Sort by ID then SIM

# ------------------------------------------------------------------------------
# Write label.data to a .csv file
	label.data.filename <- "label_simulation.csv"
	write.csv(label.data,file = label.data.filename,na = ".",quote = F,row.names = F)
