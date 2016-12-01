# in silico infliximab dosing project
# Script for simulating concentrations in the maintenance phase
# Everyone will receive 5 mg/kg dose
# ------------------------------------------------------------------------------
# Create a data frame ready for mrgsolve simulation
	# Function for creating a data frame ready for mrgsolve simulation
		label.function <- function(first.int.data) {
			# Make all predicted concentrations (IPRE) and PK parameter values after day 98 == NA
				conc.data <- first.int.data
				conc.data$IPRE[conc.data$time > max(sample.times)] <- NA
			# Assign dosing and frequency for the label simulation
				mgkg.dose <- 5	# Always 5 mg/kg
				label.int <- 56	# Always every 56 days for label simulation

			# Make a loop that pulls the patient's weight at time of next dose
			# Calculate the new dose and then administers (i.e., simulates the concentration profile)
				repeat {
					# Time of next dose
						next.dose <- max(sample.times)
					# Previous covariate values
					 	prev.WT <- conc.data$WTCOV[conc.data$time == next.dose]
						prev.ADA <- conc.data$ADA[conc.data$time == next.dose]
						prev.ALB <- conc.data$ALBCOV[conc.data$time == next.dose]
					# Calculate new dose for patient
						new.dose <- mgkg.dose*prev.WT

					# Create input data frame for simulation
						input.sim.data <- conc.data
						input.sim.data$amt[input.sim.data$time == next.dose] <- new.dose	# Add new dose to data frame at time of last sample
						input.sim.data$TIME_WT <- prev.WT
						input.sim.data$TIME_ADA <- prev.ADA
						input.sim.data$TIME_ALB <- prev.ALB
					# Re-add evid and rate columns
						input.sim.data$cmt <- 1	# Signifies which compartment the dose goes into (compartment 1 = central)
						input.sim.data$evid <- 1	# Signifies dosing event
						input.sim.data$evid[input.sim.data$amt == 0] <- 0
						input.sim.data$rate <- -2	# Signifies that infusion duration is specified in model file
						input.sim.data$rate[input.sim.data$amt == 0] <- 0
					# Flag if we want covariates to change depending on concentrations during simulation
						input.sim.data$FLAG <- time.dep
					# Simulate
						conc.data <- mod %>% mrgsim(data = input.sim.data,carry.out = c("amt","ERRPRO")) %>% as.tbl

					# Add the "next.sample" time to the list of sample.times
						next.sample <- next.dose+label.int
						sample.times <- sort(c(unique(c(sample.times,next.sample))))
					# Make all predicted concentrations (IPRE) and PK parameter values after last sample == NA
						conc.data$IPRE[conc.data$time > max(sample.times)] <- NA
					# If the last predicted concentration in the data frame (i.e., when time = 602) is NA, then continue with the loop
						if (is.na(conc.data$IPRE[conc.data$time == last.time]) == FALSE) break
				}	# Brackets closing "repeat"
			conc.data
		}
# Simulate the maintenance phase from "first.int.data1" using "label.function" for each ID (individual) in each SIM (simulation group)
	label.data <- ddply(first.int.data1, .(ID,SIM), label.function)

# ------------------------------------------------------------------------------
# Write label.data to a .csv file
	label.data.filename <- paste0("time_dep_",time.dep,"_label_simulation.csv")
	write.csv(label.data,file = label.data.filename,na = ".",quote = F,row.names = F)
