# in silico infliximab dosing project
# Script for simulating concentrations for the second, third and fourth intervals
# Everyone will receive 5 mg/kg dose
# ------------------------------------------------------------------------------
# Create a data frame ready for mrgsolve simulation
	# Function for creating a data frame ready for mrgsolve simulation
		label.function <- function(pop.data) {
			ID.number <- pop.data$ID[1]	# Individual ID
			BASE_WT <- pop.data$BASE_WT[1]	# Individual weight
			BASE_ALB <- pop.data$BASE_ALB[1]	# Individual albumin
			ETA1 <- pop.data$ETA1
			ETA2 <- pop.data$ETA2
			ETA3 <- pop.data$ETA3
			ETA4 <- pop.data$ETA4
			ERRPRO <- pop.data$ERRPRO

			input.label.data <- data.frame(
				ID = ID.number,
				time = TIME,	# Time points for simulation
				BASE_WT,	# Baseline weight
				BASE_ALB,	# Baseline albumin
				ETA1,
				ETA2,
				ETA3,
				ETA4,
				ERRPRO,
				amt = amt.init*BASE_WT,	# mg/kg dose
				evid = 1,	# Dosing event
				cmt = 1,	# Dose into the central compartment (compartment = 1)
				rate = -2	# Infusion duration is specific in the model file
			)
			# Make the amt given in the last time-point == 0
				input.label.data$amt[!c(input.label.data$time %in% TIMEi)] <- 0
				input.label.data$evid[!c(input.label.data$time %in% TIMEi)] <- 0
				input.label.data$rate[!c(input.label.data$time %in% TIMEi)] <- 0
			# Flag that this is simulation and want covariates to change depending on concentrations
				input.label.data$FLAG <- 0
		# Simulate concentration-time profiles for individuals in input.label.data
			label.data <- mod %>% mrgsim(data = input.label.data,carry.out = c("amt","ERRPRO")) %>% as.tbl
	}
# Create population data frame ready for mrgsolve simulation
	label.data <- ddply(pop.data, .(ID,SIM), label.function, .parallel = FALSE)

# ------------------------------------------------------------------------------
# Write label.data to a .csv file
	# label.data.filename <- "label_simulation.csv"
	# write.csv(label.data,file = label.data.filename,na = ".",quote = F,row.names = F)
