# Time-weighted Bayes project
# Script for running the "infliximab_bayes_simulation.R" file and then saving the output
# ------------------------------------------------------------------------------
# Create method specific results folder
	method.output.dir <- paste0(sim.output.dir,method,covariate)
	dir.create(file.path(method.output.dir),showWarnings = FALSE)
	setwd(file.path(method.output.dir))

# Source and run the Bayes scenario script
	source(paste0(work.dir,"bayes.R"))

# ------------------------------------------------------------------------------
# Save method specific results
	optimise.bayes.data.filename <- paste0(method,covariate,"_optimise_bayes_simulation.csv")
	write.csv(optimise.bayes.data,file = optimise.bayes.data.filename,na = ".",quote = F,row.names = F)

# ------------------------------------------------------------------------------
# Remove optimise.bayes.data from the workspace
	rm(optimise.bayes.data)

# ------------------------------------------------------------------------------
# Print completed message
	print(paste0("Scenario (",method," ",covariate,") has finished."))
	setwd(sim.output.dir)	# Set working directory back to global directory
