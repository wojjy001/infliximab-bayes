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
			# first.int.data <- first.int.data1
			conc.data <- first.int.data
			conc.data$IPRE[conc.data$time > max(sample.times)] <- NA
			# Previous dose mg/kg
				prev.mgkg.dose <- 5	# Initially 5 mg/kg
				prev.dose <- head(conc.data$amt,1)
				prev.int <- 56	# Initially every 56 days
				increments <- 0

		# If the last predicted concentration in the data frame (i.e., when time = 546) is NA, then continue with the loop
			repeat {
				# Time of most recent samples
					last.two.samples <- tail(sample.times,2)
					last.sample <- max(sample.times)
				# Previous DV
					prev.DV <- conc.data$DV[conc.data$time %in% last.two.samples]
					prev.DV[prev.DV < 0] <- .Machine$double.eps
				# Previous covariate values
					prev.WT <- conc.data$WTCOV[conc.data$time == last.sample]
					prev.ADA <- conc.data$ADA[conc.data$time == last.sample]
					prev.ALB <- conc.data$ALBCOV[conc.data$time == last.sample]
				# Previous dose
					prev.dose.time <- head(last.two.samples,1)
				# Calculate the new dose for the next interval based on "sample" and "dose"
					if (prev.DV[2] < 1 & prev.mgkg.dose == 5) {
						prev.mgkg.dose <- 7.5	# Now increase to 7.5 mg/kg
						new.dose <- prev.mgkg.dose*prev.WT+sum(increments)
						prev.int <- 42
					} else if (prev.DV[1] < 1 & prev.mgkg.dose == 7.5 & prev.int == 42) {
						prev.mgkg.dose <- 7.5
						new.dose <- prev.mgkg.dose*prev.WT+sum(increments)
						prev.int <- 56	# Revert back to 8-weekly dosing after the 6-week interval
					} else if (prev.DV[2] < 1 & prev.mgkg.dose == 7.5 & prev.int == 56) {
						prev.mgkg.dose <- 10
						new.dose <- prev.mgkg.dose*prev.WT+sum(increments)
						prev.int <- 42
					} else if (prev.DV[2] < 1 & prev.mgkg.dose == 10 & prev.int == 42) {
						prev.mgkg.dose <- 10
						new.dose <- prev.mgkg.dose*prev.WT+sum(increments)
						prev.int <- 28
					} else if (prev.DV[2] < 3 & prev.DV[2] >= 1) {
						increments <- c(increments,50)
						new.dose <- prev.dose+50
					} else if (prev.DV[2] > 5) {
						new.dose <- prev.dose
					} else if (prev.DV[1] > 5 & prev.DV[2] > 5) {
						increments <- c(increments,-50)
						new.dose <- prev.dose-50
					} else {
						new.dose <- prev.dose
					}
					if (new.dose > amt.max*prev.WT) {
						new.dose <- amt.max*prev.WT
						prev.mgkg.dose <- 10
					}
					prev.dose <- new.dose

				# Create input data frame for simulation
					input.sim.data <- conc.data
					input.sim.data$amt[input.sim.data$time == last.sample] <- new.dose	# Add new dose to data frame at time of last sample
					input.sim.data$TIME_WT <- prev.WT
					input.sim.data$TIME_ADA <- prev.ADA
					input.sim.data$TIME_ALB <- prev.ALB
					# Re-add evid and rate columns
						input.sim.data$cmt <- 1	# Signifies which compartment the dose goes into
						input.sim.data$evid <- 1	# Signifies dosing event
						input.sim.data$evid[input.sim.data$amt == 0] <- 0
						input.sim.data$rate <- -2	# Signifies that infusion duration is specified in model file
						input.sim.data$rate[input.sim.data$amt == 0] <- 0
					# Flag that this is simulation and want covariates to change depending on concentrations
						input.sim.data$FLAG <- time.dep
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
	clinical.data.filename <- paste0("time_dep_",time.dep,"_clinical_simulation.csv")
	write.csv(clinical.data,file = clinical.data.filename,na = ".",quote = F,row.names = F)
