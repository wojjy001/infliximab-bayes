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

# # -------------------------------------------------------------------------------
# # Parallelise jobs to increase speed
# 	library(doParallel)	# Parallel processing
# 	# Set up cores to run parallel processes, thus increasing speed
# 	# Set up a cluster of cores to run the job overall
# 		cl <- makePSOCKcluster(2)
# 		# detectCores() searches for the number of cores that the local machine has
# 	# List packages required to be sent to each core for the parallel process
# 	# The foreach package always needs to be included
# 		clusterEvalQ(cl,list(
# 			library(foreach),
# 			# source("D:/infliximab-bayes/Project/first_interval.R")	# Windows directory
# 			# source("E:/Wojciechowski/infliximab-bayes/Project/first_interval.R")	# Server directory
# 			source("/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/first_interval.R")	# Mac directory
# 		))
# 	# Register the parallel backend with the foreach package
# 		registerDoParallel(cl)

# ------------------------------------------------------------------------------
# seeds <- c(4033,96433,112209,141915,229433,299247,309449,355108,402385,416828,460403,488160,541147,604233,642803,643087,715276,725943,734678,756111,765270,968450)
# for (i in 1:length(seeds)) {
# 	seed <- seeds[i]
# Run "single-run" simulation files
# First standard interval simulation (initial dose is 5 mg/kg)
	suppressPackageStartupMessages(	# Suppress package loading messages
		suppressWarnings(	# Suppress warning messages
			source(paste0(work.dir,"first_interval1.R"))
		)
	)
# First standard interval simulation (initial dose is 10 mg/kg)
	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"first_interval2.R"))))
# Label simulation
	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"label.R"))))
# Clinical simulation where doses are adjusted based on trough concentrations (DV)
	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"clinical_TDM.R"))))
# Clinical simulation
	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"clinical.R"))))
# Run the various Bayes estimation scenarios sequentially
# Scenarios with no time-weighting
	suppressWarnings(source(paste0(work.dir,"bayes.R")))
# }
