# Time-weighted Bayes project
# Script for setting the working directory and executing other R scripts
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))
# Global directory (where R scripts are saved)
	work.dir <- "D:/infliximab-bayes/Project/"

# # -------------------------------------------------------------------------------
# # Parallelise jobs to increase speed
# 	library(doParallel)	# Parallel processing
# 	# Set up cores to run parallel processes, thus increasing speed
# 	# Set up a cluster of cores to run the job Overall
# 		cl <- makePSOCKcluster(2)
# 		# detectCores() searches for the number of cores that the local machine has
# 	# List packages required to be sent to each core for the parallel process
# 	# The foreach package always needs to be included
# 		clusterEvalQ(cl,list(
# 			library(foreach),
# 			source("D:/infliximab-bayes/Project/first_interval.R")
# 		))
# 	# Register the parallel backend with the foreach package
# 		registerDoParallel(cl)

# ------------------------------------------------------------------------------
# Run "single-run" simulation files
# First standard interval simulation
	suppressPackageStartupMessages(	# Suppress package loading messages
		suppressWarnings(	# Suppress warning messages
			source(paste0(work.dir,"first_interval.R"))
		)
	)
# Label simulation
	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"label.R"))))
# Clinical simulation
	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"clinical.R"))))

# ------------------------------------------------------------------------------
# Run the various Bayes estimation scenarios sequentially
# Scenarios with no time-weighting
	suppressWarnings(source(paste0(work.dir,"bayes.R")))
