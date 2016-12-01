# in silico infliximab dosing project
# Script for setting the working directory and executing other R scripts
# There is also a working directory for output that needs to be set in "functions.R"
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))

# Global directory (where R scripts are saved)
	work.dir <- "E:/Wojciechowski/infliximab-bayes/Project/"
# Source and compile model file
	source(paste0(work.dir,"model.R"))

# ------------------------------------------------------------------------------
# Run 12 groups of 10 simulations with both time-dependent and non-time-dependent random and covariate effect changes
	for (i in 1:12) {
		# Source universal functions file
			source(paste0(work.dir,"functions.R"))
		# Create population
			source(paste0(work.dir,"population.R"))
		# Run simulations
			try(
				for (i in 0:1) {
					time.dep <- i	# 0 = no time-dependent covariate and random effect changes
					# 1 = time-dependent covariate and random effect changes
					# Calculate ETA values for all time-points (based on time-dependence scenario)
						pop.data <- ddply(pop.data, .(SIM,ID), eta.function)
						# Write pop.data to a .csv file
							pop.data.filename <- paste0("time_dep_",time.dep,"_population_characteristics.csv")
							write.csv(pop.data,file = pop.data.filename,na = ".",quote = F,row.names = F)
					# Read in the previous pop.data
						prev.pop.data.name <- paste0("time_dep_",i,"_population_characteristics.csv")
						pop.data <- read.csv(file = prev.pop.data.name)
					# First standard interval simulation (initial dose is 5 mg/kg)
					# Common to all maintenance phase dosing strategies
						suppressWarnings(	# Suppress warning messages
							source(paste0(work.dir,"first_interval1.R"))
						)
					# Run "single-run" simulation files
					# Label simulation
						suppressWarnings(source(paste0(work.dir,"label.R")))
					# Clinical simulation where doses are adjusted based on trough concentrations (DV)
						suppressWarnings(source(paste0(work.dir,"clinical_TDM.R")))
					# Clinical protocol simulation
						suppressWarnings(source(paste0(work.dir,"clinical.R")))
					# Bayesian simulation
						suppressWarnings(source(paste0(work.dir,"bayes.R")))
				}
			)
	}
