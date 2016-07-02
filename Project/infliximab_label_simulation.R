# Time-weighted Bayes project
# Script for simulating concentrations for the second, third and fourth intervals
# Everyone will receive 5 mg/kg dose
# ------------------------------------------------------------------------------
# Source the other R scripts and execute
	source(paste0(work.dir,"infliximab_population.R"))

# ------------------------------------------------------------------------------
# Set the dose for simulating the remainder of the intervals
	amt2 <- 5	# 5 mg/kg

# Create a data frame ready for mrgsolve simulation
	# Function for creating a data frame ready for mrgsolve simulation
		population.interval2.function <- function(input.data) {
			ID.number <- input.data$ID[1]	# Individual ID
			SIM.number <- input.data$SIM[1]	# Individual simulation number
			ALB <- input.data$ALB[input.data$TIME %in% TIME]	# Individual albumin
			ADA <- input.data$ADA[input.data$TIME %in% TIME]	# Individual ADA status
			ETA1 <- input.data$ETA1[input.data$TIME %in% TIME]
			ETA2 <- input.data$ETA2[input.data$TIME %in% TIME]
			ETA3 <- input.data$ETA3[input.data$TIME %in% TIME]
			ETA4 <- input.data$ETA4[input.data$TIME %in% TIME]
			input.conc.data2 <- data.frame(
				ID = ID.number,
				SIM = SIM.number,
				time = TIME,	# Time points for simulation
				ALB,	# Albumin
				ADA,	# Anti-drug antibodies
				ETA1,
				ETA2,
				ETA3,
				ETA4,
				amt = amt2*70,	# mg/kg dose
				evid = 1,	# Dosing event
				cmt = 1,	# Dose into the central compartment (compartment = 1)
				rate = -2	# Infusion duration is specific in the model file
			)
			# Make the amt given in the last time-point == 0
			input.conc.data2$amt[!c(input.conc.data2$time %in% TIMEi)] <- 0
			input.conc.data2$evid[!c(input.conc.data2$time %in% TIMEi)] <- 0
			input.conc.data2$rate[!c(input.conc.data2$time %in% TIMEi)] <- 0
			# Return input.conc.data
				input.conc.data2
		}
	# Create population data frame ready for mrgsolve simulation
		input.conc.data2 <- ddply(pop.data, .(ID,SIM), population.interval2.function)

# Simulate concentration-time profiles for individuals in input.conc.data
	per.simulation.function2 <- function(input.data) {
		conc.data2 <- mod %>% data_set(input.data) %>% carry.out(SIM,amt) %>% mrgsim()
		conc.data2 <- as.data.frame(conc.data2)
	}
	conc.data2 <- ddply(input.conc.data2, .(SIM), per.simulation.function2)

# ------------------------------------------------------------------------------
# Test plot
	plotobj2 <- NULL
	plotobj2 <- ggplot(conc.data2)
	plotobj2 <- plotobj2 + stat_summary(aes(x = time,y = IPRE),geom = "line",fun.y = median,colour = "red")
	plotobj2 <- plotobj2 + stat_summary(aes(x = time,y = IPRE),geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "red",alpha = 0.3)
	plotobj2 <- plotobj2 + scale_y_log10("Infliximab Concentration (mg/L)\n",breaks = c(0.001,0.01,0.1,1,10,100,100),labels = c(0.001,0.01,0.1,1,10,100,100))
	plotobj2 <- plotobj2 + scale_x_continuous("\nTime (days)")
	plotobj2
