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
# Read in seeds from previous output
	project.dir <- "E:/Wojciechowski/Moved-Infliximab-Output/"
  file.list <- list.files(path = project.dir,pattern = "SUCCESS") # List of successful
	library(stringr)
	library(plyr)
	library(dplyr)
	num.list <- str_extract_all(file.list,pattern = "[0-9]")	# Find all the numbers in the filename
	numlist2 <- llply(num.list,function(num.list) paste(num.list,collapse = ""))	# Paste their numbers back together
	num.unlist <- unlist(numlist2)	# Unlist them
	num.unlist.split <- str_split(num.unlist,pattern = "[9]",n = 2)	# Where the first 9 occurs, split the string into 2
	file.data <- structure(data.frame(matrix(unlist(num.unlist.split),length(num.unlist),2,T)))
	names(file.data)[c(1,2)] <- c("nsim","seed")
	file.data$nsim <- as.numeric(levels(file.data$nsim))[file.data$nsim]
	file.data$seed <- as.numeric(levels(file.data$seed)[file.data$seed])

	for (i in 1:nrow(file.data)) {
		nsim <- file.data$nsim[i]
		seed <- file.data$seed[i]
	# Run "single-run" simulation files
	# First standard interval simulation (initial dose is 5 mg/kg)
		suppressPackageStartupMessages(	# Suppress package loading messages
			suppressWarnings(	# Suppress warning messages
				source(paste0(work.dir,"first_interval1.R"))
			)
		)
	# # First standard interval simulation (initial dose is 10 mg/kg)
	# 	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"first_interval2.R"))))
	# Label simulation
		suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"label.R"))))
	# # Clinical simulation where doses are adjusted based on trough concentrations (DV)
	# 	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"clinical_TDM.R"))))
	# # Clinical simulation
	# 	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"clinical.R"))))
	# # Run the various Bayes estimation scenarios sequentially
	# # Scenarios with no time-weighting
	# 	suppressWarnings(source(paste0(work.dir,"bayes.R")))
	}
