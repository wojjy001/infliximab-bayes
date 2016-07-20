# Time-weighted Bayes project
# Script for simulating concentrations for the second, third and fourth intervals
# Subsequent doses are dependent on measured trough concentrations
# If trough target = 3 mg/L, and measured trough is 1.5 mg/L, then next dose will be doubled
# Assuming linear kinetics - double the dose, double the trough concentration
# ------------------------------------------------------------------------------
# Simulate intervals separately
# Change doses based on "standard clinical practice" methods
	conc.data.x <- conc.data[conc.data$time < 98,]
	interval.clinical <- function(interval) {
		# Call simulated data for the previous interval
			if (interval == 2) {
				clinical.data <- conc.data	# From the first interval
				prev.TIMEXi <- TIME1i
				TIMEX <- TIME2
				TIMEXi <- TIME2i
				sample.time <- sample.time1
			} else if (interval == 3) {
				clinical.data <- clinical.data2	# From the second interval
				prev.TIMEXi <- TIME2i
				TIMEX <- TIME3
				TIMEXi <- TIME3i
				sample.time <- sample.time2
			} else {
				clinical.data <- clinical.data3	# From the third interval
				prev.TIMEXi <- TIME3i
				TIMEX <- TIME4
				TIMEXi <- TIME4i
				sample.time <- sample.time3
			}

		population.clinical <- function(input.data) {
			ID.number <- input.data$ID[1]	# Individual ID
			SIM.number <- input.data$SIM[1]	# Individual simulation number

			# clinical.data = simulated concentration data from the previous interval
				ind.clinical.data <- clinical.data[clinical.data$ID == ID.number & clinical.data$SIM == SIM.number,]
			# Pull out the sampled concentration from the individual's simulated concentration profile
				err <- ind.clinical.data$ERRPRO[ind.clinical.data$time == sample.time]	# Individual's residual error
				sample <- ind.clinical.data$IPRE[ind.clinical.data$time == sample.time]*(1+err)
				if (sample < 0) sample <- 0.0001
			# Pull out the dose that was given that resulted in that sampled concentration
				prev.dose <- ind.clinical.data$amt[ind.clinical.data$time == prev.TIMEXi[1]]

			# Calculate the new dose for the next interval based on "sample" and "dose"
				if (sample < trough.target | sample >= trough.upper) {
					new.dose <- trough.target/sample*prev.dose	# Adjust the dose if out of range
				} else {
					new.dose <- prev.dose	# Continue with previous dose if within range
				}

			# Pull the amount in the compartments at the end of the previous interval
				prev.cent <- ind.clinical.data$CENT[ind.clinical.data$time == sample.time]
				prev.peri <- ind.clinical.data$PERI[ind.clinical.data$time == sample.time]
				prev.aut <- ind.clinical.data$AUT[ind.clinical.data$time == sample.time]

			# Set up the new input data frame for mrgsolve for the next interval
				# Subset "pop.data" for the individual's data
					ind.data <- pop.data[pop.data$ID == ID.number & pop.data$SIM == SIM.number,]
				# Then call on parameter values and put into input.clinical.data
					ALB <- ind.data$ALB[ind.data$TIME %in% TIMEX]	# Individual albumin
					ADA <- ind.data$ADA[ind.data$TIME %in% TIMEX]	# Individual ADA status
					ETA1 <- ind.data$ETA1[ind.data$TIME %in% TIMEX]
					ETA2 <- ind.data$ETA2[ind.data$TIME %in% TIMEX]
					ETA3 <- ind.data$ETA3[ind.data$TIME %in% TIMEX]
					ETA4 <- ind.data$ETA4[ind.data$TIME %in% TIMEX]
					ERRPRO <- ind.data$ERRPRO[ind.data$TIME %in% TIMEX]
					input.clinical.data <- data.frame(
						ID = ID.number,
						SIM = SIM.number,
						time = TIMEX,	# Time points for simulation
						ALB,	# Albumin
						ADA,	# Anti-drug antibodies
						ETA1,
						ETA2,
						ETA3,
						ETA4,
						ERRPRO,
						amt = new.dose,	# Clinically optimised dose
						evid = 1,	# Dosing event
						cmt = 1,	# Dose into the central compartment (compartment = 1)
						rate = -2	# Infusion duration is specific in the model file
					)
				# Make the amt given in the last time-point == 0
				# Change evid and rate accordingly
					if (interval != 4) {
						input.clinical.data$amt[!c(input.clinical.data$time %in% TIMEXi) | input.clinical.data$time == max(TIMEXi)] <- 0
						input.clinical.data$evid[!c(input.clinical.data$time %in% TIMEXi) | input.clinical.data$time == max(TIMEXi)] <- 0
						input.clinical.data$rate[!c(input.clinical.data$time %in% TIMEXi) | input.clinical.data$time == max(TIMEXi)] <- 0
					} else {
						input.clinical.data$amt[!c(input.clinical.data$time %in% TIMEXi)] <- 0
						input.clinical.data$evid[!c(input.clinical.data$time %in% TIMEXi)] <- 0
						input.clinical.data$rate[!c(input.clinical.data$time %in% TIMEXi)] <- 0
					}
			# Simulate concentration time profile
				initial.compartment <- list(CENT = prev.cent,PERI = prev.peri,AUT = prev.aut)
				new.clinical.data <- mod %>% init(initial.compartment) %>% data_set(input.clinical.data) %>% carry.out(SIM,amt,ERRPRO) %>% mrgsim()
				new.clinical.data <- as.data.frame(new.clinical.data)
		}
		# Simulate concentration-time profiles for individuals in input.clinical.data
			new.clinical.data <- ddply(ID.data, .(SIM,ID), population.clinical)
	}
# Simulate the second interval
	clinical.data2 <- interval.clinical(2)
	clinical.data2.x <- clinical.data2[clinical.data2$time < 210,]
# Simulate the third interval
	clinical.data3 <- interval.clinical(3)
	clinical.data3.x <- clinical.data3[clinical.data3$time < 378,]
# Simulate the fourth interval
	clinical.data4 <- interval.clinical(4)

# Combine clinical.dataX
	clinical.data <- rbind(conc.data.x,clinical.data2.x,clinical.data3.x,clinical.data4)
	clinical.data <- clinical.data[with(clinical.data, order(clinical.data$ID,clinical.data$SIM)), ]	# Sort by ID then SIM

# ------------------------------------------------------------------------------
# Write clinical.data to a .csv file
	clinical.data.filename <- "clinical_simulation.csv"
	write.csv(clinical.data,file = clinical.data.filename,na = ".",quote = F,row.names = F)
