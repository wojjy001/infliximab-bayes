# Time-weighted Bayes project
# Script for setting the working directory and executing other R scripts
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))
# Global directory (where R scripts are saved)
	work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/"

# -------------------------------------------------------------------------------
# Create and set working directory
# Specific for the simulation
	n <- 12	# Number of seed individuals (where each seed individual has a different set of covariate values)
	nsim <- 1	# Number of simulations of the seed individuals to perform
	sim.name <- paste("SIM",nsim,"_IND",n,sep = "")	# Simulation folder's name
	sim.output.dir <- paste0(work.dir,sim.name,"/")	# Simulation directory
	dir.create(file.path(sim.output.dir),showWarnings = FALSE) # Create simulation directory
	setwd(file.path(sim.output.dir))	#Set the working directory

# -------------------------------------------------------------------------------
# Parallelise jobs to increase speed
	library(doParallel)	# Parallel processing
	# Set up cores to run parallel processes, thus increasing speed
	# Set up a cluster of cores to run the job Overall
		cl <- makePSOCKcluster(2)
		# detectCores() searches for the number of cores that the local machine has
	# List packages required to be sent to each core for the parallel process
	# The foreach package always needs to be included
		clusterEvalQ(cl,list(
			library(foreach),
			source("/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/infliximab_functions.R"),
			source("/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/infliximab_model.R"),
			source("/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/infliximab_population.R"),
			source("/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/infliximab_first_interval_simulation.R")
		))
	# Register the parallel backend with the foreach package
		registerDoParallel(cl)

# -------------------------------------------------------------------------------
# Run "single-run" simulation files
# First standard interval simulation
	suppressPackageStartupMessages(	# Suppress package loading messages
		suppressWarnings(	# Suppress warning messages
			source(paste0(work.dir,"infliximab_first_interval_simulation.R"))
		)
	)
# Label simulation
	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"infliximab_label_simulation.R"))))
# Clinical simulation
	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"infliximab_clinical_simulation.R"))))
# Optimise simulation
	suppressPackageStartupMessages(suppressWarnings(source(paste0(work.dir,"infliximab_optimise_simulation.R"))))

# -------------------------------------------------------------------------------
# Run the various Bayes estimation scenarios sequentially
# Scenarios with no time-weighting
	method <- "NTimeWeight"
	covariate <- "AllCov"	# AllCov = All covariates
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))
	covariate <- "NoADA"	# NoADA = No ADA information, assume population typical, i.e., 0
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))
	covariate <- "NoALB"	# NoALB = No albumin information, assume population typical, i.e., 4
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))
	covariate <- "NoCov"	# NoCov = No covariates, assume population typical
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))

# Scenarios using Peck, Q = 1.005
	method <- "Peck1.005"
	covariate <- "AllCov"	# AllCov = All covariates
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))
	covariate <- "NoADA"	# NoADA = No ADA information, assume population typical, i.e., 0
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))
	covariate <- "NoALB"	# NoALB = No albumin information, assume population typical, i.e., 4
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))
	covariate <- "NoCov"	# NoCov = No covariates, assume population typical
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))

# Scenarios using Peck, Q = 1.01
	method <- "Peck1.01"
	covariate <- "AllCov"	# AllCov = All covariates
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))
	covariate <- "NoADA"	# NoADA = No ADA information, assume population typical, i.e., 0
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))
	covariate <- "NoALB"	# NoALB = No albumin information, assume population typical, i.e., 4
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))
	covariate <- "NoCov"	# NoCov = No covariates, assume population typical
	suppressWarnings(source(paste0(work.dir,"infliximab_run_save.R")))

# # -------------------------------------------------------------------------------
# # Quit R once all completed
# 	q("no")
