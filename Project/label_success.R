# in silico infliximab dosing project
# Script for searching through simulation folders and labels them as SUCCESS or FAIL depending on which files exist
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))
	# system(command = "open -n -a R")

# Load package libraries
  library(stringr)

# Global directory (where R scripts are saved)
  # project.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/Moved-Infliximab-Output/" # Mac directory
  project.dir <- "E:/Wojciechowski/Moved-Infliximab-Output/"  # Server directory
  setwd(file.path(project.dir))

# Remove the "SUCCESS_" appendage to folder names if already there
# Next section will just re-add it
  success.list <- list.files(path = project.dir,pattern = "SUCCESS")
  for (i in 1:length(success.list)) {
    split.name <- str_split(success.list[i],pattern = "SUCCESS_")
    unlist.split.name <- unlist(split.name)
    file.rename(from = success.list[i],to = unlist.split.name[2])
  }

# Collate a list of the folders in the output directory
  file.list <- list.files(path = project.dir,pattern = "SIM") # List of simulation folders
# For each folder in the directory...
  for (i in 1:length(file.list)) {
    folder.dir <- paste0(project.dir,file.list[i])
    setwd(file.path(folder.dir))  # Set working directory to simulation folder directory
    if (file.exists("time_dep_0_optimise_bayes_data1.csv") == TRUE &
      file.exists("time_dep_1_optimise_bayes_data1.csv") == TRUE
    ) {
      # If the Bayes optimisation files exists, then set return to the project's directory and rename the simulation folder with "SUCCESS"
        setwd(file.path(project.dir))
        file.rename(from = file.list[i],to = paste0("SUCCESS_",file.list[i]))
    } else {
      # If the Bayes optimisation files DO NOT exist, then do nothing and just return to the project directory
        setwd(file.path(project.dir))
        # file.rename(from = test.file,to = paste0("FAIL_",test.file))
    }
  }
