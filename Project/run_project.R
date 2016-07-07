# Time-weighted Bayes project
# Script for setting the working directory and executing other R scripts
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))
# Global directory (where R scripts are saved)
	work.dir <- "D:/infliximab-bayes/Project/"

# -------------------------------------------------------------------------------
# Parallelise jobs to increase speed
	library(doParallel)	# Parallel processing
	# Set up cores to run parallel processes, thus increasing speed
	# Set up a cluster of cores to run the job Overall
		cl <- makePSOCKcluster(4)
		# detectCores() searches for the number of cores that the local machine has
	# List packages required to be sent to each core for the parallel process
	# The foreach package always needs to be included
		clusterEvalQ(cl,list(
			library(foreach),
			source("D:/infliximab-bayes/Project/first_interval.R")
		))
	# Register the parallel backend with the foreach package
		registerDoParallel(cl)

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
# # Optimise simulation
# 	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"optimise.R"))))

# ------------------------------------------------------------------------------
# Run the various Bayes estimation scenarios sequentially
# Scenarios with no time-weighting
	method <- "NTimeWeight"
	covariate <- "AllCov"	# AllCov = All covariates
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))
	covariate <- "NoADA"	# NoADA = No ADA information, assume population typical, i.e., 0
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))
	covariate <- "NoALB"	# NoALB = No albumin information, assume population typical, i.e., 4
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))
	covariate <- "NoCov"	# NoCov = No covariates, assume population typical
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))

# Scenarios using Peck, Q = 1.005
	method <- "Peck1.005"
	covariate <- "AllCov"	# AllCov = All covariates
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))
	covariate <- "NoADA"	# NoADA = No ADA information, assume population typical, i.e., 0
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))
	covariate <- "NoALB"	# NoALB = No albumin information, assume population typical, i.e., 4
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))
	covariate <- "NoCov"	# NoCov = No covariates, assume population typical
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))

# Scenarios using Peck, Q = 1.01
	method <- "Peck1.01"
	covariate <- "AllCov"	# AllCov = All covariates
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))
	covariate <- "NoADA"	# NoADA = No ADA information, assume population typical, i.e., 0
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))
	covariate <- "NoALB"	# NoALB = No albumin information, assume population typical, i.e., 4
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))
	covariate <- "NoCov"	# NoCov = No covariates, assume population typical
	suppressWarnings(source(paste0(work.dir,"run_bayes_save.R")))

# -----------------------------------------------------------------------------
# Quit R once all completed
	q("no")
