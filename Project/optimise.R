# Time-weighted Bayes project
# Script for simulating concentrations for the second, third and fourth intervals
# Subsequent doses are dependent on measured trough concentrations
# Doses are optimised using maximum likelihood estimation
# ------------------------------------------------------------------------------
# Call simulated data for the first interval
	work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/"
	source(paste0(work.dir,"infliximab_first_interval_simulation.R"))

	optimise.data1 <- conc.data

# Optimise doses using maximum likelihood estimation
	interval.optimise <- function(prev.TIMEXi,TIMEX,TIMEXi,sample.time) {

		if (max(prev.TIMEXi) == 98) {
			optimise.data <- optimise.data1
		} else if (max(prev.TIMEXi) == 210) {
			optimise.data <- optimise.data2
		} else {
			optimise.data <- optimise.data3
		}

		population.optimise <- function(ID.data) {
			ID.number <- ID.data$ID[1]	# Individual ID
			SIM.number <- ID.data$SIM[1]	# Individual simulation number

			# optimise.data = simulated concentration data from the previous interval
				ind.optimise.data <- optimise.data[optimise.data$ID == ID.number & optimise.data$SIM == SIM.number,]
			# Pull out the sampled concentration from the individual's simulated concentration profile
				prev.IPRE <- ind.optimise.data$IPRE[ind.optimise.data$time == sample.time]
			# Pull ou thte dose from the previous interval
				prev.dose <- ind.optimise.data$amt[ind.optimise.data$time == prev.TIMEXi[1]]

			# Pull the amount in the compartments at the end of the previous interval
				prev.cent <- ind.optimise.data$CENT[ind.optimise.data$time == sample.time]
				prev.peri <- ind.optimise.data$PERI[ind.optimise.data$time == sample.time]
				prev.aut <- ind.optimise.data$AUT[ind.optimise.data$time == sample.time]

			# Set up the new input data frame for mrgsolve for the next interval
				# Subset "pop.data" for the individual's data
					ind.data <- pop.data[pop.data$ID == ID.number & pop.data$SIM == SIM.number & pop.data$TIME %in% TIMEX,]
				# Then call on parameter values and put into input.optimise.data
					ALB <- ind.data$ALB	# Individual albumin
					ADA <- ind.data$ADA	# Individual ADA status
					ETA1 <- ind.data$ETA1
					ETA2 <- ind.data$ETA2
					ETA3 <- ind.data$ETA3
					ETA4 <- ind.data$ETA4
					ERRPRO <- ind.data$ERRPRO
					input.optimise.data <- data.frame(
						ID = ID.number,
						SIM = SIM.number,
						time = TIMEX,	# Time points for simulation
						ALB,	# Albumin
						ADA,	# Anti-drug antibodies
						ETA1,
						ETA2,
						ETA3,
						ETA4,
						ERRPRO,
						amt = prev.dose,	# optimised dose
						evid = 1,	# Dosing event
						cmt = 1,	# Dose into the central compartment (compartment = 1)
						rate = -2	# Infusion duration is specific in the model file
					)
					# Make the amt given in the last time-point == 0
					# Change evid and rate accordingly
						if (max(TIMEXi) != 490) {
							input.optimise.data$amt[!c(input.optimise.data$time %in% TIMEXi) | input.optimise.data$time == max(TIMEXi)] <- 0
							input.optimise.data$evid[!c(input.optimise.data$time %in% TIMEXi) | input.optimise.data$time == max(TIMEXi)] <- 0
							input.optimise.data$rate[!c(input.optimise.data$time %in% TIMEXi) | input.optimise.data$time == max(TIMEXi)] <- 0
						}
						# If looking at the fourth interval, still give the last time-point a dose
							if (max(TIMEXi) == 490) {
								input.optimise.data$amt[!c(input.optimise.data$time %in% TIMEXi)] <- 0
								input.optimise.data$evid[!c(input.optimise.data$time %in% TIMEXi)] <- 0
								input.optimise.data$rate[!c(input.optimise.data$time %in% TIMEXi)] <- 0
							}

			# Modify model code ready for simulation
				initial.compartment <- list(CENT = prev.cent,PERI = prev.peri,AUT = prev.aut)
				mod <- mod %>% init(initial.compartment) %>% carry.out(SIM,amt,ERRPRO)
			# Reduce input data frame to only infusion and trough times
			# This will be the data frame that will be used for dose optimisation - as few time-points as possible
			# optim function is a bit of a bottleneck
				trough.times <- TIMEXi+56	# Infusion time + 56 days
				optimise.times <- unique(sort(c(TIMEXi,trough.times)))
				input.optimise.dose.data <- input.optimise.data[input.optimise.data$time %in% optimise.times,]

				# Initial parameter estimates
					# Initial dose is what the "clinical" dose would have been - somewhere in the ballpark
					# Pull out the sampled concentration from the individual's simulated concentration profile
						err <- ind.optimise.data$ERRPRO[ind.optimise.data$time == sample.time]	# Individual's residual error
						sample <- ind.optimise.data$IPRE[ind.optimise.data$time == sample.time]*exp(err)
					# Pull out the dose that was given that resulted in that sampled concentration
						prev.dose <- ind.optimise.data$amt[ind.optimise.data$time == prev.TIMEXi[1]]

					# Calculate the new dose for the next interval based on "sample" and "dose"
						if (sample < trough.target | sample >= trough.upper) {
							initial.dose <- trough.target/sample*prev.dose	# Adjust the dose if out of range
						} else {
							initial.dose <- prev.dose	# Continue with previous dose if within range
						}

					initial.err <- 0.001	# Initial estimate for the error
					initial.par <- c(initial.dose,initial.dose,initial.dose,initial.err)
					par <- initial.par

				# Dose parameter limits
					lower.dose.limit <- initial.dose*0.5
					upper.dose.limit <- initial.dose*2

				# Find the doses that maximum the likelihood of trough concentrations being the target
					optimise.dose <- function(par) {
						# Add fitted parameters to the input data frame
							input.optimise.dose.data$amt[input.optimise.dose.data$time == TIMEXi[1]] <- par[1]
							input.optimise.dose.data$amt[input.optimise.dose.data$time == TIMEXi[2]] <- par[2]
							input.optimise.dose.data$amt[input.optimise.dose.data$time == TIMEXi[3]] <- par[3]
						# Simulate concentration-time profiles with fitted doses
							new.optimise.data <- mod %>% data_set(input.optimise.dose.data) %>% mrgsim()
							new.optimise.data <- as.data.frame(new.optimise.data)
						# Pull out the predicted trough concentrations with the fitted doses for the interval
							yhat <- new.optimise.data$IPRE
							res <- dnorm(trough.target,yhat,yhat*par[4],log = T)	# Minimise the error between target trough (3 mg/L) and predicted trough concentrations
						# Objective function value and minimise the value
							objective <- -1*sum(res)
							objective
					}
					optimised.doses <- optim(par,optimise.dose,hessian = FALSE,method = "L-BFGS-B",lower = c(lower.dose.limit,lower.dose.limit,lower.dose.limit,0.0001),upper = c(upper.dose.limit,upper.dose.limit,upper.dose.limit,1))
					optimised.doses
				# Input optimised doses ready for simulation
					input.optimise.data$amt[input.optimise.data$time == TIMEXi[1]] <- optimised.doses$par[1]	# First dose
					input.optimise.data$amt[input.optimise.data$time == TIMEXi[2]] <- optimised.doses$par[2]	# Second dose
					input.optimise.data$amt[input.optimise.data$time == TIMEXi[3]] <- optimised.doses$par[3]	# Third dose
			# Simulate
				new.optimise.data <- mod %>% data_set(input.optimise.data) %>% mrgsim()
				new.optimise.data <- as.data.frame(new.optimise.data)
		}
		# Simulate concentration-time profiles for individuals in ID.data
			new.optimise.data <- ddply(ID.data, .(SIM,ID), population.optimise, .parallel = TRUE, .progress = "text")
			new.optimise.data
	}

# Simulate intervals separately
# Simulate the second interval
	optimise.data2 <- interval.optimise(TIME1i,TIME2,TIME2i,sample.time1)
# Simulate the third interval
	optimise.data3 <- interval.optimise(TIME2i,TIME3,TIME3i,sample.time2)
# Simulate the fourth interval
	optimise.data4 <- interval.optimise(TIME3i,TIME4,TIME4i,sample.time3)

# Combine optimise.dataX
	optimise.data <- rbind(optimise.data1,optimise.data2,optimise.data3,optimise.data4)
	optimise.data <- optimise.data[with(optimise.data, order(optimise.data$ID,optimise.data$SIM)), ]	# Sort by ID then SIM

# ------------------------------------------------------------------------------
# Write optimise.data to a .csv file
	optimise.data.filename <- "optimise_simulation.csv"
	write.csv(optimise.data,file = optimise.data.filename,na = ".",quote = F,row.names = F)

# # ------------------------------------------------------------------------------
# # Test plot
# 	plotobj4 <- NULL
# 	plotobj4 <- ggplot(optimise.data)
# 	plotobj4 <- plotobj4 + stat_summary(aes(x = time,y = IPRE),geom = "line",fun.y = median,colour = "red")
# 	plotobj4 <- plotobj4 + stat_summary(aes(x = time,y = IPRE),geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "red",alpha = 0.3)
# 	plotobj4 <- plotobj4 + geom_hline(aes(yintercept = trough.target),linetype = "dashed")
# 	plotobj4 <- plotobj4 + geom_hline(aes(yintercept = trough.upper),linetype = "dashed")
# 	plotobj4 <- plotobj4 + scale_y_log10("Infliximab Concentration (mg/L)\n",breaks = c(0.001,0.01,0.1,1,10,100,100),labels = c(0.001,0.01,0.1,1,10,100,100))
# 	plotobj4 <- plotobj4 + scale_x_continuous("\nTime (days)")
# 	plotobj4
