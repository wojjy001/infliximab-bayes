# Time-weighted Bayes project
# Script for setting the working directory and executing other R scripts
#-------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))
	
# Set working directory
	work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/"
	setwd(work.dir)

#-------------------------------------------------------------------------------
# Source the other R scripts and execute
	# infliximab_population.R
	source(paste0(work.dir,"infliximab_population.R"))
