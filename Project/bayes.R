# Time-weighted Bayes project
# Script for estimating individual parameter values and using them to guide dose optimisation
# Subsequent doses are dependent on measured trough concentrations
# Doses are optimised using maximum likelihood estimation
# ------------------------------------------------------------------------------
# Optimise doses using maximum likelihood estimation
	conc.data.x <- conc.data[conc.data$time < 98,]
	interval.bayes.optimise <- function(interval) {
		# Call simulated data for the previous interval
			if (interval == 2) {
				optimise.bayes.data <- conc.data	# From the first interval
				TIMEX <- unique(c(TIME1,TIME2))	# Time sequence for the first x intervals and the next one
				prev.TIMEX <- TIME1	# Time sequence for previous x intervals
				next.TIMEX <- TIME2	# Next time interval's sequence
				TIMEXi <- unique(c(TIME1i,TIME2i))	# Infusion times for the previous x intervals and the next one
				prev.TIMEXi <- TIME1i	# Infusion times from the previous x intervals
				next.TIMEXi <- TIME2i	# Infusion times for the next interval
				sample.times <- sample.time1	# Sample times from the previous x intervals
				sample.time <- sample.time1	# Most recent sample time
			} else if (interval == 3) {
				optimise.bayes.data <- rbind(conc.data.x,optimise.bayes.data2)	# From the first two intervals
				TIMEX <- unique(c(TIME1,TIME2,TIME3))
				prev.TIMEX <- unique(c(TIME1,TIME2))
				next.TIMEX <- TIME3
				TIMEXi <- unique(c(TIME1i,TIME2i,TIME3i))
				prev.TIMEXi <- unique(c(TIME1i,TIME2i))
				next.TIMEXi <- TIME3i
				sample.times <- c(sample.time1,sample.time2)
				sample.time <- sample.time2
			} else {
				optimise.bayes.data <- rbind(conc.data.x,optimise.bayes.data2.x,optimise.bayes.data3)	# From the first three intervals
				TIMEX <- unique(c(TIME1,TIME2,TIME3,TIME4))
				prev.TIMEX <- unique(c(TIME1,TIME2,TIME3))
				next.TIMEX <- TIME4
				TIMEXi <- unique(c(TIME1i,TIME2i,TIME3i,TIME4i))
				prev.TIMEXi <- unique(c(TIME1i,TIME2i,TIME3i))
				next.TIMEXi <- TIME4i
				sample.times <- c(sample.time1,sample.time2,sample.time3)
				sample.time <- sample.time3
			}
		# Call in other objects/data frames required
		# Need to be placed here for parallelisation to work
			covariate.scenario <- covariate	# Level of covariate information for scenario
			method.scenario <- method	# Time-weighting method
			pop.data.import <- pop.data	# Data frame of population's characteristics

		###########
		##_BAYES_##
		###########
		# Estimate individual parameter values using:
			# 1. Trough sample from the end of the previous interval
			# 2. Covariate values measured at beginning and end of previous interval
			# 3. Known doses that were administered during the previous interval

		# Create an input data frame for Bayesian estimation AND THEN simulation given the results of  Bayesian estimation
			input.bayes <- function(ID.data) {
				ID.number <- ID.data$ID[1]	# Individual ID
				SIM.number <- ID.data$SIM[1]	# Individual simulation number

				# optimise.bayes.data = simulated concentration data from the previous interval
					prev.bayes.data <- optimise.bayes.data[optimise.bayes.data$ID == ID.number & optimise.bayes.data$SIM == SIM.number,]
				# Pull out the dose from the previous interval
					prev.dose <- prev.bayes.data$amt[prev.bayes.data$amt != 0]

				# Pull covariate information from the beginning and end of the previous interval
					if (covariate.scenario == "AllCov" | covariate.scenario == "NoADA") {
						ALB1 <- prev.bayes.data$ALB[prev.bayes.data$time == 0]
						ALB2 <- prev.bayes.data$ALB[prev.bayes.data$time == sample.time1]
						ALB3 <- prev.bayes.data$ALB[prev.bayes.data$time == sample.time2]
						ALB4 <- prev.bayes.data$ALB[prev.bayes.data$time == sample.time3]
						if (interval == 2) {
							TIMEalb <- c(0,sample.time1)
							RATEalb <- c(ALB1,ALB2)
						}
						if (interval == 3) {
							TIMEalb <- c(0,sample.time1,sample.time2)
							RATEalb <- c(ALB1,ALB2,ALB3)
						}
						if (interval == 4) {
							TIMEalb <- c(0,sample.time1,sample.time2,sample.time3)
							RATEalb <- c(ALB1,ALB2,ALB3,ALB4)
						}
						extrap.alb <- approxfun(TIMEalb,RATEalb,method = "linear")
						prev.ALB <- extrap.alb(prev.TIMEX)
					}
					if (covariate.scenario == "NoALB" | covariate.scenario == "NoCov") {
						prev.ALB <- 4
					}
					if (covariate.scenario == "AllCov" | covariate.scenario == "NoALB") {
						ADA1 <- prev.bayes.data$ADA[prev.bayes.data$time == sample.time1]
						ADA2 <- prev.bayes.data$ADA[prev.bayes.data$time == sample.time2]
						ADA3 <- prev.bayes.data$ADA[prev.bayes.data$time == sample.time3]
						if (interval == 2) {
							RATEada <- c(ADA1,ADA1)
						}
						if (interval == 3) {
							RATEada <- c(ADA1,ADA2,ADA2)
						}
						if (interval == 4) {
							RATEada <- c(ADA1,ADA2,ADA3,ADA3)
						}
						TIMEada <- c(0,sample.times)
						extrap.ada <- approxfun(TIMEada,RATEada,method = "constant")
						prev.ADA <- extrap.ada(prev.TIMEX)
					}
					if (covariate.scenario == "NoADA" | covariate.scenario == "NoCov") {
						prev.ADA <- 0
					}

				# Set up the new input data frame for Bayes estimation
					input.bayes.data <- data.frame(
						ID = ID.number,
						SIM = SIM.number,
						time = prev.TIMEX,	# Time points for simulation
						ALB = prev.ALB,	# Albumin
						ADA = prev.ADA,	# Anti-drug antibodies
						ETA1 = 0,
						ETA2 = 0,
						ETA3 = 0,
						ETA4 = 0,
						ERRPRO = 0,
						amt = 350,	# Automatically writes in the first dose, will replace other doses for other intervals
						evid = 1,	# Dosing event
						cmt = 1,	# Dose into the central compartment (compartment = 1)
						rate = -2	# Infusion duration is specific in the model file
					)
				# Input doses according to the interval they were actually administered in
					if (interval > 2) input.bayes.data$amt[input.bayes.data$time %in% TIME2i] <- c(prev.dose[4],prev.dose[5],0)
					if (interval > 3) input.bayes.data$amt[input.bayes.data$time %in% TIME3i] <- c(prev.dose[6],prev.dose[7],prev.dose[8],0)
				# Make the amt given in the last time-point == 0 - only sampling here, not dosing straight away
				# Change evid and rate accordingly
					input.bayes.data$amt[!c(input.bayes.data$time %in% prev.TIMEXi) | input.bayes.data$time >= max(prev.TIMEXi)] <- 0
					input.bayes.data$evid[!c(input.bayes.data$time %in% prev.TIMEXi) | input.bayes.data$time >= max(prev.TIMEXi)] <- 0
					input.bayes.data$rate[!c(input.bayes.data$time %in% prev.TIMEXi) | input.bayes.data$time >= max(prev.TIMEXi)] <- 0
				input.bayes.data
			}
			input.bayes.data <- ddply(ID.data, .(SIM,ID), input.bayes, .parallel = TRUE, .progress = "text")
		# Write results to .csv
			filename1 <- paste0(method.scenario,covariate.scenario,"_interval",interval,"_input_bayes.csv")
			write.csv(input.bayes.data,file = filename1,na = ".",quote = FALSE,row.names = FALSE)

		# Estimate individual parameter values using "input.bayes.data"
		# Note: It does not contain the previous samples
		# Should have the previous doses
			bayes.eta <- function(ID.data) {
				ID.number <- ID.data$ID[1]	# Individual ID
				SIM.number <- ID.data$SIM[1]	# Individual simulation number

				# optimise.bayes.data = simulated concentration data from the previous interval
					prev.bayes.data <- optimise.bayes.data[optimise.bayes.data$ID == ID.number & optimise.bayes.data$SIM == SIM.number,]
				# Pull out the sampled concentration from the individual's simulated concentration profile
					prev.err <- prev.bayes.data$ERRPRO[prev.bayes.data$time %in% sample.times]
					prev.DV <- prev.bayes.data$IPRE[prev.bayes.data$time %in% sample.times]*exp(prev.err)

				# Subset input.bayes.data for ID and SIM
				# Only collect the previous infusion times - speed up estimation process
					input.bayes.data <- input.bayes.data[input.bayes.data$ID == ID.number & input.bayes.data$SIM == SIM.number & input.bayes.data$time %in% prev.TIMEXi,]

				# Initial parameter estimates
					initial.bayes.par <- exp(c(0,0,0,0))
					par <- initial.bayes.par
				# Bayesian estimation
					bayes.estimate <- function(par) {
						# Describe parameters to be optimised
							ETA1fit <- log(par[1])	# In the exponential domain
							ETA2fit <- log(par[2])
							ETA3fit <- log(par[3])
							ETA4fit <- log(par[4])
						# Add fitted parameters to the input data frame
							input.bayes.data$ETA1 <- ETA1fit
							input.bayes.data$ETA2 <- ETA2fit
							input.bayes.data$ETA3 <- ETA3fit
							input.bayes.data$ETA4 <- ETA4fit
						# Simulate concentration-time profiles with fitted parameters
							new.bayes.data <- mod %>% data_set(input.bayes.data) %>% mrgsim()
							new.bayes.data <- as.data.frame(new.bayes.data)
						# Pull out the predicted trough concentrations with the fitted doses for the interval
							yhat <- new.bayes.data$IPRE[new.bayes.data$time %in% sample.times]
							# Posterior log-likelihood
							# Error model: Y = IPRE*exp(ERRPRO), log(Y) = log(IPRE) + ERRPRO
								TIMET <- max(new.bayes.data$time) - new.bayes.data$time	# Time since last observation
								if (method.scenario == "NTimeWeight") loglikpost.sd = ERRPRO	# No time-weighting
								if (method.scenario == "Peck1.005") loglikpost.sd = ERRPRO*1.005^TIMET  # Peck method, Q = 1.005
								if (method.scenario == "Peck1.01")	loglikpost.sd = ERRPRO*1.01^TIMET	# Peck method, Q = 1.01
								loglikpost <- dnorm(log(prev.DV),mean = log(yhat),sd = loglikpost.sd,log = T)
							# Prior log-likelihood
								ETA <- c(ETA1fit,ETA2fit,ETA3fit,ETA4fit)
								ETABSV <- c(PPVCL,PPVV1,PPVQ,PPVV2)
								loglikprior <- dnorm(ETA,mean = 0,sd = ETABSV,log = T)
						# Objective function value and minimise the value
							objective <- -1*sum(loglikpost,loglikprior)
					}
				# Run bayes.estimate function
					bayes.result <- optim(par,bayes.estimate,hessian = FALSE,method = "L-BFGS-B",lower = c(0.001,0.001,0.001,0.001),upper = c(Inf,Inf,Inf,Inf),control = list(parscale = par,factr = 1e7))

				# Convert new ETA values (estimated in the exp() domain)
					new.ETA1 <- log(bayes.result$par[1])
					new.ETA2 <- log(bayes.result$par[2])
					new.ETA3 <- log(bayes.result$par[3])
					new.ETA4 <- log(bayes.result$par[4])

				# Output data frame of individually fitted ETAs
					bayes.eta.data <- data.frame(
						ID = ID.number,
						SIM = SIM.number,
						interval,
						ETA1 = new.ETA1,
						ETA2 = new.ETA2,
						ETA3 = new.ETA3,
						ETA4 = new.ETA4
					)
			}
			bayes.eta.data <- ddply(ID.data, .(SIM,ID), bayes.eta, .parallel = TRUE, .progress = "text")
		# Write results to .csv
			filename2 <- paste0(method.scenario,covariate.scenario,"_interval",interval,"_bayes_eta.csv")
			write.csv(bayes.eta.data,file = filename2,na = ".",quote = FALSE,row.names = FALSE)

		# Simulate the Bayesian estimated profile given results in "bayes.eta.data" and using "input.bayes.data" as a template
			bayes.sim <- function(ID.data) {
				ID.number <- ID.data$ID[1]	# Individual ID
				SIM.number <- ID.data$SIM[1]	# Individual simulation number

				# Subset input.bayes.data for ID and SIM
					input.bayes.data <- input.bayes.data[input.bayes.data$ID == ID.number & input.bayes.data$SIM == SIM.number,]
				# Subset bayes.eta.data for ID and SIM
					bayes.eta.data <- bayes.eta.data[bayes.eta.data$ID == ID.number & bayes.eta.data$SIM == SIM.number,]
				# Input estimated individual parameter values and covariate values for simulation of previous interval
					input.bayes.data$ETA1 <- bayes.eta.data$ETA1[1]	# ETA for clearance
					input.bayes.data$ETA2 <- bayes.eta.data$ETA2[1]	# ETA for V1
					input.bayes.data$ETA3 <- bayes.eta.data$ETA3[1]	# ETA for Q
					input.bayes.data$ETA4 <- bayes.eta.data$ETA4[1]	# ETA for V2
				# Simulate previous interval according to the individual parameter estimates and level of covariate information
					bayes.sim.data <- mod %>% data_set(input.bayes.data) %>% mrgsim()
					bayes.sim.data <- as.data.frame(bayes.sim.data)
			}
			bayes.sim.data <- ddply(ID.data, .(SIM,ID), bayes.sim, .parallel = TRUE, .progress = "text")
		# Write results to .csv
			filename3 <- paste0(method.scenario,covariate.scenario,"_interval",interval,"_bayes_sim.csv")
			write.csv(bayes.sim.data,file = filename3,na = ".",quote = FALSE,row.names = FALSE)

		##############
		##_OPTIMISE_##
		##############
		# Optimise dose for the individual using:
			# 1. Individual BAYES predicted concentration (compartment amounts)
			# 2. Carried forward covariate information from the end of the previous interval
		# OR:
			# 1. Individual IPRE predicted concentration
			# 2. Actual covariate information

		# Optimise dose for the next interval given the estimated Bayes parameters and carried forward covariate values
		# Begins from end of bayes.sim.data
			bayes.optim <- function(ID.data) {
				ID.number <- ID.data$ID[1]	# Individual ID
				SIM.number <- ID.data$SIM[1]	# Individual simulation number

				# Subset bayes.sim.data for ID and SIM
					bayes.sim.data <- bayes.sim.data[bayes.sim.data$ID == ID.number & bayes.sim.data$SIM == SIM.number,]
				# Subset bayes.eta.data for ID and SIM
					bayes.eta.data <- bayes.eta.data[bayes.eta.data$ID == ID.number & bayes.eta.data$SIM == SIM.number,]

				# Call on last covariate values from the previous interval
				# This values are carried forward for the next interval
					last.ALB <- bayes.sim.data$ALB[bayes.sim.data$time == sample.time]	# Albumin
					last.ADA <- bayes.sim.data$ADA[bayes.sim.data$time == sample.time]	# Anti-drug antibodies
				# Call on previously estimated Bayes parameter values
					bayes.ETA1 <- bayes.eta.data$ETA1[1]
					bayes.ETA2 <- bayes.eta.data$ETA2[1]
					bayes.ETA3 <- bayes.eta.data$ETA3[1]
					bayes.ETA4 <- bayes.eta.data$ETA4[1]

				# Set up a data frame to be input for dose optimisation based on Bayes parameters
					input.optimise.data <- data.frame(
						ID = ID.number,
						SIM = SIM.number,
						time = next.TIMEXi,	# Time-points for Simulation (set for every 14 days)
						ALB = last.ALB[1],	# Carry forward last albumin concentration
						ADA = last.ADA[1],	# Carry forward last ADA status
						ETA1 = bayes.ETA1,
						ETA2 = bayes.ETA2,
						ETA3 = bayes.ETA3,
						ETA4 = bayes.ETA4,
						ERRPRO = 0,
						amt = 0,	# Dose to be optimised
						evid = 1,	# Dosing event
						cmt = 1,	# Dose into the central compartment (compartment = 1)
						rate = -2	# Infusion duration is specific in the model file
					)
				# Make the amt given in the last time-point == 0
				# Change evid and rate accordingly
					if (interval != 4) {
						input.optimise.data$amt[!c(input.optimise.data$time %in% next.TIMEXi) | input.optimise.data$time == max(next.TIMEXi)] <- 0
						input.optimise.data$evid[!c(input.optimise.data$time %in% next.TIMEXi) | input.optimise.data$time == max(next.TIMEXi)] <- 0
						input.optimise.data$rate[!c(input.optimise.data$time %in% next.TIMEXi) | input.optimise.data$time == max(next.TIMEXi)] <- 0
					}
				# If looking at the fourth interval, still give the last time-point a dose
					if (interval == 4) {
						input.optimise.data$amt[!c(input.optimise.data$time %in% next.TIMEXi)] <- 0
						input.optimise.data$evid[!c(input.optimise.data$time %in% next.TIMEXi)] <- 0
						input.optimise.data$rate[!c(input.optimise.data$time %in% next.TIMEXi)] <- 0
					}

				# Pull the amount in the compartments at the end of the previous interval
					prev.bayes.cent <- bayes.sim.data$CENT[bayes.sim.data$time == sample.time]
					prev.bayes.peri <- bayes.sim.data$PERI[bayes.sim.data$time == sample.time]
					prev.bayes.aut <- bayes.sim.data$AUT[bayes.sim.data$time == sample.time]
				# Modify model code ready for simulation
					initial.bayes.compartment <- list(CENT = prev.bayes.cent,PERI = prev.bayes.peri,AUT = prev.bayes.aut)
					optim.mod1 <- mod %>% init(initial.bayes.compartment) %>% carry.out(SIM,amt,ERRPRO)

				# Subset optimise.bayes.data for ID and SIM
					optimise.bayes.data <- optimise.bayes.data[optimise.bayes.data$ID == ID.number & optimise.bayes.data$SIM == SIM.number,]
				# Initial parameter estimates
					# Initial dose is what the "clinical" dose would have been - somewhere in the ballpark
					# Based on last sample and last dose
						last.err <- optimise.bayes.data$ERRPRO[optimise.bayes.data$time == sample.time]
						last.DV <- optimise.bayes.data$IPRE[optimise.bayes.data$time == sample.time]*exp(last.err)
						last.dose <- optimise.bayes.data$amt[optimise.bayes.data$time == max(prev.TIMEXi[prev.TIMEXi != max(prev.TIMEXi)])]
					# Calculate the new dose for the next interval based on "sample" and "dose"
						if (last.DV < trough.target | last.DV >= trough.upper) {
							initial.dose <- trough.target/last.DV*last.dose	# Adjust the dose if out of range
						} else {
							initial.dose <- last.dose	# Continue with previous dose if within range
						}
					# Limits of parameters
						if (interval == 2) {
							initial.par <- c(initial.dose,initial.dose,0.001)
							lower.limit <- c(initial.dose*0.5,initial.dose*0.5,0.0001)
							upper.limit <- c(initial.dose*2,initial.dose*2,10)
						} else {
							initial.par <- c(initial.dose,initial.dose,initial.dose,0.001)
							lower.limit <- c(initial.dose*0.5,initial.dose*0.5,initial.dose*0.5,0.0001)
							upper.limit <- c(initial.dose*2,initial.dose*2,initial.dose*2,10)
						}
						par <- initial.par

				# Reduce input data frame to only infusion and trough times
				# This will be the data frame that will be used for dose optimisation - as few time-points as possible
				# optim function is a bit of a bottleneck
					trough.times <- next.TIMEXi+56	# Infusion time + 56 days
					optimise.times <- unique(sort(c(next.TIMEXi,trough.times)))
					input.optimise.data <- input.optimise.data[input.optimise.data$time %in% optimise.times,]

				# Find the doses that maximum the likelihood of trough concentrations being the target
					optimise.dose <- function(par) {
						# Add fitted parameters to the input data frame
							input.optimise.data$amt[input.optimise.data$time == next.TIMEXi[1]] <- par[1]
							input.optimise.data$amt[input.optimise.data$time == next.TIMEXi[2]] <- par[2]
							if (interval == 2) {
								err <- par[3]
							} else {
								input.optimise.data$amt[input.optimise.data$time == next.TIMEXi[3]] <- par[3]
								err <- par[4]
							}
						# Simulate concentration-time profiles with fitted doses
							new.optimise.data <- optim.mod1 %>% data_set(input.optimise.data) %>% mrgsim()
							new.optimise.data <- as.data.frame(new.optimise.data)
						# Pull out the predicted trough concentrations with the fitted doses for the interval
							yhat <- new.optimise.data$IPRE
							res <- dnorm(trough.target,yhat,yhat*err,log = T)	# Minimise the error between target trough (3 mg/L) and predicted trough concentrations
						# Objective function value and minimise the value
							objective <- -1*sum(res)
							objective
					}
					optimised.doses <- optim(par,optimise.dose,hessian = FALSE,method = "L-BFGS-B",lower = lower.limit,upper = upper.limit)

				# Create a data frame of results
					bayes.optim.data <- data.frame(
						ID = ID.number,
						SIM = SIM.number,
						interval,
						DOSE1 = optimised.doses$par[1],
						DOSE2 = optimised.doses$par[2],
						DOSE3 = NA
					)
					if (interval != 2) {	# Other intervals contain three doses
						bayes.optim.data$DOSE3 = optimised.doses$par[3]
					}
					bayes.optim.data
			}
			bayes.optim.data <- ddply(ID.data, .(SIM,ID), bayes.optim, .parallel = TRUE, .progress = "text")
		# Write to .csv
			filename4 <- paste0(method.scenario,covariate.scenario,"_interval",interval,"_bayes_optim.csv")
			write.csv(bayes.optim.data,file = filename4,na = ".",quote = FALSE,row.names = FALSE)

		# Optimise dose for the next interval given the true individual parameters and true covariate values
		# Begins from end of optimise.bayes.data
			sim.optim <- function(ID.data) {
				ID.number <- ID.data$ID[1]	# Individual ID
				SIM.number <- ID.data$SIM[1]	# Individual simulation number

				# Subset optimise.bayes.data for ID and SIM
					optimise.bayes.data <- optimise.bayes.data[optimise.bayes.data$ID == ID.number & optimise.bayes.data$SIM == SIM.number,]
				# Call on previous amounts in the compartments, doses and samples
					# Sampled concentration from the individual's simulated concentration profile
						last.err <- optimise.bayes.data$ERRPRO[optimise.bayes.data$time == sample.time]	# Individual's residual error
						sample <- optimise.bayes.data$IPRE[optimise.bayes.data$time == sample.time]*exp(last.err)
					# First dose from the previous interval
						last.dose <- optimise.bayes.data$amt[optimise.bayes.data$time == prev.TIMEXi[1]]
					# Amount in the compartments at the end of the previous interval
						last.cent <- optimise.bayes.data$CENT[optimise.bayes.data$time == sample.time]
						last.peri <- optimise.bayes.data$PERI[optimise.bayes.data$time == sample.time]
						last.aut <- optimise.bayes.data$AUT[optimise.bayes.data$time == sample.time]
					# Modify model code ready for simulation
						initial.compartment <- list(CENT = last.cent,PERI = last.peri,AUT = last.aut)
						optim.mod2 <- mod %>% init(initial.compartment) %>% carry.out(SIM,amt,ERRPRO)

				# Call on individual's true covariate and random effect values from pop.data.import
					# Subset "pop.data.import" for ID and SIM
						ind.data <- pop.data.import[pop.data.import$ID == ID.number & pop.data.import$SIM == SIM.number & pop.data.import$TIME %in% TIMEX,]
					# Covariate values
						ALB <- ind.data$ALB	# Individual albumin
						ADA <- ind.data$ADA	# Individual ADA status
					# Random effect values
						ETA1 <- ind.data$ETA1
						ETA2 <- ind.data$ETA2
						ETA3 <- ind.data$ETA3
						ETA4 <- ind.data$ETA4
						ERRPRO <- ind.data$ERRPRO

				# Set up a data frame to be input for dose optimisation based on true parameters
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
						amt = last.dose,	# Optimised dose
						evid = 1,	# Dosing event
						cmt = 1,	# Dose into the central compartment (compartment = 1)
						rate = -2	# Infusion duration is specific in the model file
					)
				# Make the amt given in the last time-point == 0
				# Change evid and rate accordingly
					if (interval != 4) {
						input.optimise.data$amt[!c(input.optimise.data$time %in% TIMEXi) | input.optimise.data$time == max(TIMEXi)] <- 0
						input.optimise.data$evid[!c(input.optimise.data$time %in% TIMEXi) | input.optimise.data$time == max(TIMEXi)] <- 0
						input.optimise.data$rate[!c(input.optimise.data$time %in% TIMEXi) | input.optimise.data$time == max(TIMEXi)] <- 0
					}
				# If looking at the fourth interval, still give the last time-point a dose
					if (interval == 4) {
						input.optimise.data$amt[!c(input.optimise.data$time %in% TIMEXi)] <- 0
						input.optimise.data$evid[!c(input.optimise.data$time %in% TIMEXi)] <- 0
						input.optimise.data$rate[!c(input.optimise.data$time %in% TIMEXi)] <- 0
					}

					# Initial parameter estimates
						# Initial dose is what the "clinical" dose would have been - somewhere in the ballpark
						# Calculate the new dose for the next interval based on "sample" and "dose"
							if (sample < trough.target | sample >= trough.upper) {
								initial.dose <- trough.target/sample*last.dose	# Adjust the dose if out of range
							} else {
								initial.dose <- last.dose	# Continue with previous dose if within range
							}
						# Limits of parameters
							if (interval == 2) {
								initial.par <- c(initial.dose,initial.dose,0.001)
								lower.limit <- c(initial.dose*0.5,initial.dose*0.5,0.0001)
								upper.limit <- c(initial.dose*2,initial.dose*2,10)
							} else {
								initial.par <- c(initial.dose,initial.dose,initial.dose,0.001)
								lower.limit <- c(initial.dose*0.5,initial.dose*0.5,initial.dose*0.5,0.0001)
								upper.limit <- c(initial.dose*2,initial.dose*2,initial.dose*2,10)
							}
							par <- initial.par

					# Reduce input data frame to only infusion and trough times
					# This will be the data frame that will be used for dose optimisation - as few time-points as possible
						trough.times <- next.TIMEXi+56	# Infusion time + 56 days
						optimise.times <- unique(sort(c(next.TIMEXi,trough.times)))
						input.optimise.data <- input.optimise.data[input.optimise.data$time %in% optimise.times,]

					# Find the doses that maximum the likelihood of trough concentrations being the target
						optimise.dose <- function(par) {
							# Add fitted parameters to the input data frame
								input.optimise.data$amt[input.optimise.data$time == next.TIMEXi[1]] <- par[1]
								input.optimise.data$amt[input.optimise.data$time == next.TIMEXi[2]] <- par[2]
								if (interval == 2) {
									err <- par[3]
								} else {
									input.optimise.data$amt[input.optimise.data$time == next.TIMEXi[3]] <- par[3]
									err <- par[4]
								}
							# Simulate concentration-time profiles with fitted doses
								new.optimise.data <- optim.mod2 %>% data_set(input.optimise.data) %>% mrgsim()
								new.optimise.data <- as.data.frame(new.optimise.data)
							# Pull out the predicted trough concentrations with the fitted doses for the interval
								yhat <- new.optimise.data$IPRE
								res <- dnorm(trough.target,yhat,yhat*err,log = T)	# Minimise the error between target trough (3 mg/L) and predicted trough concentrations
							# Objective function value and minimise the value
								objective <- -1*sum(res)
								objective
						}
						optimised.doses <- optim(par,optimise.dose,hessian = FALSE,method = "L-BFGS-B",lower = lower.limit,upper = upper.limit)

				# Create a data frame of results
					sim.optim.data <- data.frame(
						ID = ID.number,
						SIM = SIM.number,
						interval,
						DOSE1 = optimised.doses$par[1],
						DOSE2 = optimised.doses$par[2],
						DOSE3 = NA
					)
					if (interval != 2) {
						sim.optim.data$DOSE3 = optimised.doses$par[3]
					}
					sim.optim.data
			}
			sim.optim.data <- ddply(ID.data, .(SIM,ID), sim.optim, .parallel = TRUE, .progress = "text")
		# Write to .csv
			filename5 <- paste0(method.scenario,covariate.scenario,"_interval",interval,"_sim_optim.csv")
			write.csv(sim.optim.data,file = filename5,na = ".",quote = FALSE,row.names = FALSE)

		# Simulate the real concentration profile for the individual given the bayes optimised doses for the  next interval
		# Begins from end of optimise.bayes.data
			real.sim <- function(ID.data) {
				ID.number <- ID.data$ID[1]	# Individual ID
				SIM.number <- ID.data$SIM[1]	# Individual simulation number

				# Subset optimise.bayes.data for ID and SIM
					optimise.bayes.data <- optimise.bayes.data[optimise.bayes.data$ID == ID.number & optimise.bayes.data$SIM == SIM.number,]
					# Amount in the compartments at the end of the previous interval
						last.cent <- optimise.bayes.data$CENT[optimise.bayes.data$time == sample.time]
						last.peri <- optimise.bayes.data$PERI[optimise.bayes.data$time == sample.time]
						last.aut <- optimise.bayes.data$AUT[optimise.bayes.data$time == sample.time]
					# Modify model code ready for simulation
						initial.compartment <- list(CENT = last.cent,PERI = last.peri,AUT = last.aut)
						optim.mod3 <- mod %>% init(initial.compartment) %>% carry.out(SIM,amt,ERRPRO)

				# Call on individual's true covariate and random effect values from pop.data.import
					# Subset "pop.data.import" for ID and SIM
						ind.data <- pop.data.import[pop.data.import$ID == ID.number & pop.data.import$SIM == SIM.number & pop.data.import$TIME %in% next.TIMEX,]
					# Covariate values
						ALB <- ind.data$ALB	# Individual albumin
						ADA <- ind.data$ADA	# Individual ADA status
					# Random effect values
						ETA1 <- ind.data$ETA1
						ETA2 <- ind.data$ETA2
						ETA3 <- ind.data$ETA3
						ETA4 <- ind.data$ETA4
						ERRPRO <- ind.data$ERRPRO

				# Set up a data frame to be input for dose optimisation based on true parameters
					input.sim.data <- data.frame(
						ID = ID.number,
						SIM = SIM.number,
						time = next.TIMEX,	# Time points for simulation
						ALB,	# Albumin
						ADA,	# Anti-drug antibodies
						ETA1,
						ETA2,
						ETA3,
						ETA4,
						ERRPRO,
						amt = 350,	# Wil be replaced with the actual doses
						evid = 1,	# Dosing event
						cmt = 1,	# Dose into the central compartment (compartment = 1)
						rate = -2	# Infusion duration is specific in the model file
					)
				# Make the amt given in the last time-point == 0
				# Change evid and rate accordingly
					if (interval != 4) {
						input.sim.data$amt[!c(input.sim.data$time %in% TIMEXi) | input.sim.data$time == max(TIMEXi)] <- 0
						input.sim.data$evid[!c(input.sim.data$time %in% TIMEXi) | input.sim.data$time == max(TIMEXi)] <- 0
						input.sim.data$rate[!c(input.sim.data$time %in% TIMEXi) | input.sim.data$time == max(TIMEXi)] <- 0
					}
				# If looking at the fourth interval, still give the last time-point a dose
					if (interval == 4) {
						input.sim.data$amt[!c(input.sim.data$time %in% TIMEXi)] <- 0
						input.sim.data$evid[!c(input.sim.data$time %in% TIMEXi)] <- 0
						input.sim.data$rate[!c(input.sim.data$time %in% TIMEXi)] <- 0
					}

				# Subset bayes.optim.data for SIM and ID
					bayes.optim.data <- bayes.optim.data[bayes.optim.data$ID == ID.number & bayes.optim.data$SIM == SIM.number & bayes.optim.data$interval == interval,]
					# Input optimised doses ready for simulation
						input.sim.data$amt[input.sim.data$time == next.TIMEXi[1]] <- bayes.optim.data$DOSE1	# First dose
						input.sim.data$amt[input.sim.data$time == next.TIMEXi[2]] <- bayes.optim.data$DOSE2	# Second dose
						if (interval != 2) input.sim.data$amt[input.sim.data$time == next.TIMEXi[3]] <- bayes.optim.data$DOSE3	# Third dose

				# Simulate
					real.sim.data <- optim.mod3 %>% data_set(input.sim.data) %>% mrgsim()
					real.sim.data <- as.data.frame(real.sim.data)
			}
			real.sim.data <- ddply(ID.data, .(SIM,ID), real.sim, .parallel = TRUE, .progress = "text")
	}

# Simulate intervals separately
# Simulate the second interval
	optimise.bayes.data2 <- interval.bayes.optimise(2)
	optimise.bayes.data2.x <- optimise.bayes.data2[optimise.bayes.data2$time < 210,]
# Simulate the third interval
	optimise.bayes.data3 <- interval.bayes.optimise(3)
	optimise.bayes.data3.x <- optimise.bayes.data3[optimise.bayes.data3$time < 378,]
# Simulate the fourth interval
	optimise.bayes.data4 <- interval.bayes.optimise(4)

# Combine bayes.dataX
	optimise.bayes.data <- rbind(conc.data.x,optimise.bayes.data2.x,optimise.bayes.data3.x,optimise.bayes.data4)
	optimise.bayes.data <- optimise.bayes.data[with(optimise.bayes.data, order(optimise.bayes.data$ID,optimise.bayes.data$SIM)), ]	# Sort by ID then SIM

# # ------------------------------------------------------------------------------
# # Test plot
# 	plotobj5 <- NULL
# 	plotobj5 <- ggplot()
# 	plotobj5 <- plotobj5 + stat_summary(aes(x = time,y = IPRE),data = optimise.bayes.data,geom = "line",fun.y = median,colour = "red")
# 	plotobj5 <- plotobj5 + stat_summary(aes(x = time,y = IPRE),data = optimise.bayes.data,geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "red",alpha = 0.3)
# 	# plotobj5 <- plotobj5 + stat_summary(aes(x = time,y = IPRE),data = optimise.data,geom = "line",fun.y = median,colour = "blue")
# 	# plotobj5 <- plotobj5 + stat_summary(aes(x = time,y = IPRE),data = optimise.data,geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "blue",alpha = 0.3)
# 	plotobj5 <- plotobj5 + geom_hline(aes(yintercept = trough.target),linetype = "dashed")
# 	plotobj5 <- plotobj5 + geom_hline(aes(yintercept = trough.upper),linetype = "dashed")
# 	plotobj5 <- plotobj5 + scale_y_log10("Infliximab Concentration (mg/L)\n",breaks = c(0.001,0.01,0.1,1,10,100,100),labels = c(0.001,0.01,0.1,1,10,100,100))
# 	plotobj5 <- plotobj5 + scale_x_continuous("\nTime (days)")
# 	plotobj5
