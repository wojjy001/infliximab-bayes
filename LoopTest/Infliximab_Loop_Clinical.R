# Script that wraps an individual patient into a simulating-dose adjustment loop
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list=ls(all=TRUE))

# Load package libraries
	library(plyr) # Split and rearrange data, ddply function
	library(dplyr)	# Split and rearrange data
	library(mrgsolve)	# Metrum Research Group differential equation solver
	library(ggplot2)	# Plotting
	library(grid)	# Plotting

# Set working directory
	work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/LoopTest/"
	setwd(work.dir)

# ------------------------------------------------------------------------------
# Source the infliximab first_interval file
	source(paste0(work.dir,"first_interval.R"))
	# Output object is called "first.int.data"
	# Sample a random individual from the simulated data
		sample.ID <- sample(1:20,1)
		# sample.ID <- 20
		conc.data <- first.int.data[first.int.data$SIM == 1 & first.int.data$ID == sample.ID,] # First three doses for the individual have been simulated

# ------------------------------------------------------------------------------
# Set up a loop that will sample the individual's concentration, estimate empirical Bayes parameters, optimise their dose and administer until time = 546 days
# Define the last time-point to be simulated
	last.time <- 546	# days
# After the initiation phase, the first sample will be collected at day 98
	sample.times <- c(0,98)	# days
# Initial dosing interval for the maintenance phase
	dose.int <- 56	# days
# Make all predicted concentrations (IPRE) and PK parameter values after sample.time1 == NA
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
			conc.data <- mod %>% carry.out(amt,SIM,ERRPRO) %>% data_set(input.sim.data) %>% mrgsim()
			conc.data <- as.data.frame(conc.data)
		# Add the "next.sample" time to the list of sample.times
			next.sample <- last.sample+dose.int
			sample.times <- sort(c(unique(c(sample.times,next.sample))))

		# Make all predicted concentrations (IPRE) and PK parameter values after sample.time1 == NA
			conc.data$IPRE[conc.data$time > max(sample.times)] <- NA

		# Plot individual's concentrations as they are being administered and optimised
			scale.log10.labels <- c(0.01,0.1,1,10,100,1000)
			plotobj1 <- NULL
			plotobj1 <- ggplot(conc.data[conc.data$time <= last.time,])
			plotobj1 <- plotobj1 + geom_line(aes(x = time,y = IPRE),colour = "red")
			plotobj1 <- plotobj1 + geom_hline(aes(yintercept = trough.target),linetype = "dashed")
			plotobj1 <- plotobj1 + geom_hline(aes(yintercept = trough.upper),linetype = "dashed")
			plotobj1 <- plotobj1 + scale_y_log10("Infliximab Concentration (mg/L)\n",breaks = scale.log10.labels,labels = scale.log10.labels,lim = c(0.01,10000))
			plotobj1 <- plotobj1 + scale_x_continuous("\nTime (days)",breaks = c(0,14,42,98,154,210,266,322,378,434,490,546))
			print(plotobj1)

		if (is.na(conc.data$IPRE[conc.data$time == last.time]) == FALSE) break

	}	# Brackets closing "repeat"

# ------------------------------------------------------------------------------
# Numerical summary of individual
# Infusion times and amounts
	print(
		data.frame(
			ID = sample.ID,
			times = conc.data$time[conc.data$amt != 0],
			amt = conc.data$amt[conc.data$amt != 0],
			int = c(0,diff(conc.data$time[conc.data$amt != 0]))
		)
	)
#	Cumulative time under target trough by time = 546
 	print(conc.data$TBT[conc.data$time == last.time])
# Print cumulative time under target trough by time = 98
	# See the effects of lame initiation phase dosing
		print(conc.data$TBT[conc.data$time == 98])
# Print cumulative area under target trough at time = 546
	print(conc.data$AUT[conc.data$time == last.time])
