# Time-weighted Bayes project
# Script for simulating concentrations for the second, third and fourth intervals
# Subsequent doses are dependent on measured trough concentrations
# If trough target = 3 mg/L, and measured trough is 1.5 mg/L, then next dose will be doubled
# Assuming linear kinetics - double the dose, double the trough concentration
# ------------------------------------------------------------------------------
# Source the other R scripts and execute
	source(paste0(work.dir,"infliximab_population.R"))

# ------------------------------------------------------------------------------
