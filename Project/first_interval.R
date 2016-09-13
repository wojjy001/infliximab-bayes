# in silico infliximab dosing project
# Script for simulating concentrations for the first interval
# Everyone will receive 5 mg/kg doses for the first interval in all scenarios
# ------------------------------------------------------------------------------
# Source the other R scripts and execute
	# work.dir <- "D:/infliximab-bayes/Project/"	# Windows directory
	# work.dir <- "E:/Wojciechowski/infliximab-bayes/Project/"	# Server directory
	work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/"	# Mac directory
	source(paste0(work.dir,"population.R"))
	source(paste0(work.dir,"model.R"))

# ------------------------------------------------------------------------------
# Create a data frame ready for mrgsolve simulation
	# Function for creating a data frame ready for mrgsolve simulation
		first.int.function <- function(pop.data) {
			ID.number <- pop.data$ID[1]	# Individual ID
			BASE_WT <- pop.data$BASE_WT[1]	# Individual weight
			BASE_ALB <- pop.data$BASE_ALB[1]	# Individual albumin
			ETA1 <- pop.data$ETA1
			ETA2 <- pop.data$ETA2
			ETA3 <- pop.data$ETA3
			ETA4 <- pop.data$ETA4
			ERRPRO <- pop.data$ERRPRO

			input.first.int.data <- data.frame(
				ID = ID.number,
				time = TIME,	# Time points for simulation
				BASE_WT,	# Baseline weight
				BASE_ALB,	# Baseline albumin
				ETA1,
				ETA2,
				ETA3,
				ETA4,
				ERRPRO,
				amt = amt.init*BASE_WT,	# mg/kg dose
				evid = 1,	# Dosing event
				cmt = 1,	# Dose into the central compartment (compartment = 1)
				rate = -2	# Infusion duration is specific in the model file
			)
			# Make the amt given in the last time-point == 0
				input.first.int.data$amt[!c(input.first.int.data$time %in% TIME1i)] <- 0
				input.first.int.data$evid[!c(input.first.int.data$time %in% TIME1i)] <- 0
				input.first.int.data$rate[!c(input.first.int.data$time %in% TIME1i)] <- 0
			# Flag that this is simulation and want covariates to change depending on concentrations
				input.first.int.data$FLAG <- 0
			# Simulate concentration-time profiles for individuals in input.conc.data
				first.int.data <- mod %>% mrgsim(data = input.first.int.data,carry.out = c("amt","ERRPRO")) %>% as.tbl
		}
# Create population data frame ready for mrgsolve simulation
	first.int.data <- ddply(pop.data, .(ID,SIM), first.int.function, .parallel = FALSE)
