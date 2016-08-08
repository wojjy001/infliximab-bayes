# Time-weighted Bayes project
# Script for simulating concentrations for the second, third and fourth intervals
# Subsequent doses are dependent on measured trough concentrations
# If trough target = 3 mg/L, and measured trough is 1.5 mg/L, then next dose will be doubled
# Assuming linear kinetics - double the dose, double the trough concentration
# ------------------------------------------------------------------------------
# Set up a loop that will sample the individual's concentration optimise their dose and administer until time = 546 days
	clinical.function <- function(ID.data) {
		SIM.number <- ID.data$SIM[1]
		ID.number <- ID.data$ID[1]
		# Make all predicted concentrations (IPRE) and PK parameter values after sample.time1 == NA
			conc.data <- first.int.data[first.int.data$SIM == SIM.number & first.int.data$ID == ID.number,]
			conc.data$IPRE[conc.data$time > max(sample.times)] <- NA

		# If the last predicted concentration in the data frame (i.e., when time = 546) is NA, then continue with the loop
			repeat {
				# Time of most recent sample
					last.sample <- max(sample.times)
				# Previous DV
					prev.DV <- conc.data$DV[conc.data$time == last.sample]
					prev.DV[prev.DV < 0] <- 0.001
				# Previous dose
					prev.dose.time <- head(tail(sample.times,2),1)
					prev.dose <- conc.data$amt[conc.data$time == prev.dose.time]
				# Calculate the new dose for the next interval based on "sample" and "dose"
					if (prev.DV < trough.target | prev.DV >= trough.upper) {
						new.dose <- trough.target/prev.DV*prev.dose	# Adjust the dose if out of range
					} else {
						new.dose <- prev.dose	# Continue with previous dose if within range
					}
				# Cap "new.dose" to 50 mg/kg (i.e., 3500 mg)
					if (new.dose > 3500) new.dose <- 3500

				# Create input data frame for simulation
					input.sim.data <- conc.data
					input.sim.data$amt[input.sim.data$time == last.sample] <- new.dose	# Add new dose to data frame at time of last sample
					# Re-add evid and rate columns
						input.sim.data$cmt <- 1	# Signifies which compartment the dose goes into
						input.sim.data$evid <- 1	# Signifies dosing event
						input.sim.data$evid[input.sim.data$amt == 0] <- 0
						input.sim.data$rate <- -2	# Signifies that infusion duration is specified in model file
						input.sim.data$rate[input.sim.data$amt == 0] <- 0
				# Simulate
					conc.data <- mod %>% mrgsim(data = input.sim.data,carry.out = c("amt","ERRPRO")) %>% as.tbl
				# Add the "next.sample" time to the list of sample.times
					next.sample <- last.sample+dose.int
					sample.times <- sort(c(unique(c(sample.times,next.sample))))
				# Make all predicted concentrations (IPRE) and PK parameter values after sample.time1 == NA
					conc.data$IPRE[conc.data$time > max(sample.times)] <- NA
				# If the last predicted concentration in the data frame (i.e., when time = 546) is NA, then continue with the loop
					if (is.na(conc.data$IPRE[conc.data$time == last.time]) == FALSE) break
			}	# Brackets closing "repeat"
		conc.data
	}	# Brackets closing "clinical.function"

	clinical.data <- ddply(ID.data, .(SIM,ID), clinical.function, .parallel = TRUE)

# ------------------------------------------------------------------------------
# Write clinical.data to a .csv file
	clinical.data.filename <- "clinical_simulation.csv"
	write.csv(clinical.data,file = clinical.data.filename,na = ".",quote = F,row.names = F)
