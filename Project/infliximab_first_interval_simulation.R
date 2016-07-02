# Time-weighted Bayes project
# Script for simulating concentrations for the first interval
# Everyone will receive 5 mg/kg doses for the first interval in all scenarios
# ------------------------------------------------------------------------------
# Source the other R scripts and execute
	source(paste0(work.dir,"infliximab_population.R"))

# ------------------------------------------------------------------------------
# Set the dose for simulating the first intervals
	amt1 <- 5	# 5 mg/kg

# Create a data frame ready for mrgsolve simulation
	# Function for creating a data frame ready for mrgsolve simulation
		population.interval1.function <- function(input.data) {
			ID.number <- input.data$ID[1]	# Individual ID
			SIM.number <- input.data$SIM[1]	# Individual simulation number
			ALB <- input.data$ALB[input.data$TIME %in% TIME1]	# Individual albumin
			ADA <- input.data$ADA[input.data$TIME %in% TIME1]	# Individual ADA status
			ETA1 <- input.data$ETA1[input.data$TIME %in% TIME1]
			ETA2 <- input.data$ETA2[input.data$TIME %in% TIME1]
			ETA3 <- input.data$ETA3[input.data$TIME %in% TIME1]
			ETA4 <- input.data$ETA4[input.data$TIME %in% TIME1]
			input.conc.data1 <- data.frame(
				ID = ID.number,
				SIM = SIM.number,
				time = TIME1,	# Time points for simulation
				ALB,	# Albumin
				ADA,	# Anti-drug antibodies
				ETA1,
				ETA2,
				ETA3,
				ETA4,
				amt = amt1*70,	# mg/kg dose
				evid = 1,	# Dosing event
				cmt = 1,	# Dose into the central compartment (compartment = 1)
				rate = -2	# Infusion duration is specific in the model file
			)
			# Make the amt given in the last time-point == 0
				input.conc.data1$amt[!c(input.conc.data1$time %in% TIME1i) | input.conc.data1$time == max(TIME1i)] <- 0
				input.conc.data1$evid[!c(input.conc.data1$time %in% TIME1i) | input.conc.data1$time == max(TIME1i)] <- 0
				input.conc.data1$rate[!c(input.conc.data1$time %in% TIME1i) | input.conc.data1$time == max(TIME1i)] <- 0
			# Return input.conc.data
				input.conc.data1
		}
	# Create population data frame ready for mrgsolve simulation
		input.conc.data1 <- ddply(pop.data, .(ID,SIM), population.interval1.function)

# Simulate concentration-time profiles for individuals in input.conc.data
	per.simulation.function1 <- function(input.data) {
		conc.data1 <- mod %>% data_set(input.data) %>% carry.out(SIM,amt) %>% mrgsim()
		conc.data1 <- as.data.frame(conc.data1)
	}
	conc.data1 <- ddply(input.conc.data1, .(SIM), per.simulation.function1)

# ------------------------------------------------------------------------------
# Test plot
	plotobj1 <- NULL
	plotobj1 <- ggplot(conc.data1)
	plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPRE),geom = "line",fun.y = median,colour = "red")
	plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPRE),geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "red",alpha = 0.3)
	plotobj1 <- plotobj1 + scale_y_log10("Infliximab Concentration (mg/L)\n",breaks = c(0.001,0.01,0.1,1,10,100,100),labels = c(0.001,0.01,0.1,1,10,100,100))
	plotobj1 <- plotobj1 + scale_x_continuous("\nTime (days)")
	plotobj1
