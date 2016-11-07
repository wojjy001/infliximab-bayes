# in silico infliximab dosing project
# Script for setting the working directory and executing other R scripts
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))
	# system(command = "open -n -a R")

# Global directory (where R scripts are saved)
	# work.dir <- "D:/infliximab-bayes/Project/"	# Windows directory
	work.dir <- "E:/Wojciechowski/infliximab-bayes/Project/"	# Server directory
	# work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/"	# Mac directory
# Source and compile model file
	source(paste0(work.dir,"model.R"))

# ------------------------------------------------------------------------------
# Read in seeds from previous output
	# project.dir <- "E:/Wojciechowski/Moved-Infliximab-Output/"
  # file.list <- list.files(path = project.dir,pattern = "SUCCESS") # List of successful
	# library(stringr)
	# library(plyr)
	# library(dplyr)
	# num.list <- str_extract_all(file.list,pattern = "[0-9]")	# Find all the numbers in the filename
	# numlist2 <- llply(num.list,function(num.list) paste(num.list,collapse = ""))	# Paste their numbers back together
	# num.unlist <- unlist(numlist2)	# Unlist them
	# num.unlist.split <- str_split(num.unlist,pattern = "[9]",n = 2)	# Where the first 9 occurs, split the string into 2
	# file.data <- structure(data.frame(matrix(unlist(num.unlist.split),length(num.unlist),2,T)))
	# names(file.data)[c(1,2)] <- c("nsim","seed")
	# file.data$nsim <- as.numeric(levels(file.data$nsim))[file.data$nsim]
	# file.data$seed <- as.numeric(levels(file.data$seed)[file.data$seed])
	#
for (i in 1:15) {
	# Source universal functions file
		source(paste0(work.dir,"functions.R"))
		# nsim <- file.data$nsim[i]
		# seed <- file.data$seed[i]
		print(paste0("seed ",seed," nsim ",nsim))
		# Create population
			source(paste0(work.dir,"population.R"))
	try(
		for (i in 0:1) {
			time.dep <- i	# 0 = no time-dependent covariate and random effect changes
			# 1 = time-dependent covariate and random effect changes
			# Calculate ETA values for all time-points (based on time-dependence scenario)
				pop.data <- ddply(pop.data, .(SIM,ID), eta.function)
				# Write pop.data to a .csv file
					pop.data.filename <- paste0("time_dep_",time.dep,"_population_characteristics.csv")
					write.csv(pop.data,file = pop.data.filename,na = ".",quote = F,row.names = F)
			# First standard interval simulation (initial dose is 5 mg/kg)
				suppressWarnings(	# Suppress warning messages
					source(paste0(work.dir,"first_interval1.R"))
				)
			# Run "single-run" simulation files
			# Label simulation
				suppressWarnings(source(paste0(work.dir,"label.R")))
			# Clinical simulation where doses are adjusted based on trough concentrations (DV)
				suppressWarnings(source(paste0(work.dir,"clinical_TDM.R")))
			# Clinical simulation
				suppressWarnings(source(paste0(work.dir,"clinical.R")))
			# Bayesian simulation
				suppressWarnings(source(paste0(work.dir,"bayes.R")))
		}
	)
}
