# in silico infliximab dosing project
# Script for simulating concentrations for the second, third and fourth intervals
# Subsequent doses are dependent on measured trough concentrations
# If trough target = 3 mg/L, and measured trough is 1.5 mg/L, then next dose will be doubled
# Assuming linear kinetics - double the dose, double the trough concentration
# ------------------------------------------------------------------------------
# Set up a loop that will sample the individual's concentration optimise their dose and administer until time = 546 days
	# first.int.data <- first.int.data[first.int.data$ID == 1 & first.int.data$SIM == 1,]
	clinical.function <- function(first.int.data) {
		# Make all predicted concentrations (IPRE) and PK parameter values after sample.time1 == NA
			conc.data <- first.int.data
			conc.data$IPRE[conc.data$time > max(sample.times)] <- NA
			# Previous dose mg/kg
				prev.mgkg.dose <- 5	# Initially 5 mg/kg
				prev.int <- 56	# Initially every 56 days

		# If the last predicted concentration in the data frame (i.e., when time = 546) is NA, then continue with the loop
			repeat {
				# Time of most recent sample
					last.sample <- max(sample.times)
				# Previous DV
					prev.DV <- conc.data$DV[conc.data$time == last.sample]
					prev.DV[prev.DV < 0] <- 1e-8
				# Previous weight
					prev.WT <- conc.data$WTCOV[conc.data$time == last.sample]
				# Previous dose
					prev.dose.time <- head(tail(sample.times,2),1)
				# Calculate the new dose for the next interval based on "sample" and "dose"
					if (prev.DV < trough.target) {
						if (prev.mgkg.dose == 5) {
							prev.mgkg.dose <- 7.5	# Now increase to 7.5 mg/kg
							new.dose <- prev.mgkg.dose*prev.WT
							prev.int <- 56
						} else if (prev.mgkg.dose == 7.5) {
							prev.mgkg.dose <- 10	# Now increase to 10 mg/kg
							new.dose <- prev.mgkg.dose*prev.WT
							prev.int <- 56
						} else {
							prev.mgkg.dose <- 10
							new.dose <- prev.mgkg.dose*prev.WT
							prev.int <- 42	# Now increase frequency to every 42 days instead of every 56 days
						}
					} else {
						new.dose <- prev.mgkg.dose*prev.WT	# Continue with previous dose if within range
					}

				# Create input data frame for simulation
					input.sim.data <- conc.data
					input.sim.data$amt[input.sim.data$time == last.sample] <- new.dose	# Add new dose to data frame at time of last sample
					# Re-add evid and rate columns
						input.sim.data$cmt <- 1	# Signifies which compartment the dose goes into
						input.sim.data$evid <- 1	# Signifies dosing event
						input.sim.data$evid[input.sim.data$amt == 0] <- 0
						input.sim.data$rate <- -2	# Signifies that infusion duration is specified in model file
						input.sim.data$rate[input.sim.data$amt == 0] <- 0
					# Flag that this is simulation and want covariates to change depending on concentrations
						input.sim.data$FLAG <- 0
				# Simulate
					conc.data <- mod %>% mrgsim(data = input.sim.data,carry.out = c("amt","ERRPRO")) %>% as.tbl
				# Add the "next.sample" time to the list of sample.times
					next.sample <- last.sample+prev.int
					sample.times <- sort(c(unique(c(sample.times,next.sample))))
				# Make all predicted concentrations (IPRE) and PK parameter values after sample.time1 == NA
					conc.data$IPRE[conc.data$time > max(sample.times)] <- NA
				# If the last predicted concentration in the data frame (i.e., when time = 546) is NA, then continue with the loop
					if (is.na(conc.data$IPRE[conc.data$time == last.time]) == FALSE) break
			}	# Brackets closing "repeat"
		conc.data
	}	# Brackets closing "clinical.function"

	clinical.data <- ddply(first.int.data1, .(SIM,ID), clinical.function, .parallel = FALSE)

# ------------------------------------------------------------------------------
# Write clinical.data to a .csv file
	clinical.data.filename <- "clinical_simulation.csv"
	write.csv(clinical.data,file = clinical.data.filename,na = ".",quote = F,row.names = F)
