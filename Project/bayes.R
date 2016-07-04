# Time-weighted Bayes project
# Script for estimating individual parameter values and using them to guide dose optimisation
# Subsequent doses are dependent on measured trough concentrations
# Doses are optimised using maximum likelihood estimation
# ------------------------------------------------------------------------------
# Optimise doses using maximum likelihood estimation
	interval.bayes <- function(prev.TIMEXi,TIMEX,TIMEXi,sample.time) {
		# Call simulated data for the previous interval
			if (max(prev.TIMEXi) == 98) {
				optimise.bayes.data <- conc.data	# From the first interval
			} else if (max(prev.TIMEXi) == 210) {
				optimise.bayes.data <- optimise.bayes.data2	# From the second interval
			} else {
				optimise.bayes.data <- optimise.bayes.data3	# From the third interval
			}
		# Call in other objects/data frames required
		# Need to be placed here for parallelisation to work
			covariate.scenario <- covariate	# Level of covariate information for scenario
			method.scenario <- method	# Time-weighting method
			pop.data.import <- pop.data	# Data frame of population's characteristics

		population.bayes <- function(ID.data) {
			ID.number <- ID.data$ID[1]	# Individual ID
			SIM.number <- ID.data$SIM[1]	# Individual simulation number

			###########
			##_BAYES_##
			###########
			# Estimate individual parameter values using:
				# 1. Trough sample from the end of the previous interval
				# 2. Covariate values measured at beginning and end of previous interval
				# 3. Known doses that were administered during the previous interval

			# optimise.bayes.data = simulated concentration data from the previous interval
				ind.bayes.data <- optimise.bayes.data[optimise.bayes.data$ID == ID.number & optimise.bayes.data$SIM == SIM.number,]
			# Pull out the sampled concentration from the individual's simulated concentration profile
				prev.err <- ind.bayes.data$ERRPRO[ind.bayes.data$time == sample.time]
				prev.DV <- ind.bayes.data$IPRE[ind.bayes.data$time == sample.time]*exp(prev.err)
			# Pull out the dose from the previous interval
				prev.dose <- ind.bayes.data$amt[ind.bayes.data$time == prev.TIMEXi[1]]

			# Pull the amount in the compartments at the beginning of the previous interval
				prev.init.cent <- ind.bayes.data$CENT[ind.bayes.data$time == prev.TIMEXi[1]]
				prev.init.peri <- ind.bayes.data$PERI[ind.bayes.data$time == prev.TIMEXi[1]]
				prev.init.aut <- ind.bayes.data$AUT[ind.bayes.data$time == prev.TIMEXi[1]]

			# Pull out previous ETA values - guide optimisation algorithm with decent initial estimates
				prev.ETA1 <- 0
				prev.ETA2 <- 0
				prev.ETA3 <- 0
				prev.ETA4 <- 0
				# If looking not looking at the first interval, make ETAs be the previously estimated ETA value
					if (prev.TIMEXi[1] != 0) {
						prev.ETA1 <- ind.bayes.data$ETA1[ind.bayes.data$time == prev.TIMEXi[1]]
						prev.ETA2 <- ind.bayes.data$ETA2[ind.bayes.data$time == prev.TIMEXi[1]]
						prev.ETA3 <- ind.bayes.data$ETA3[ind.bayes.data$time == prev.TIMEXi[1]]
						prev.ETA4 <- ind.bayes.data$ETA4[ind.bayes.data$time == prev.TIMEXi[1]]
					}

			# Pull covariate information from at the beginning and end of the previous interval
				if (covariate.scenario == "AllCov") {	# All covariate information is used
					BASE_ALB <- ind.bayes.data$ALB[ind.bayes.data$time == prev.TIMEXi[1]]
					FINAL_ALB <- ind.bayes.data$ALB[ind.bayes.data$time == tail(prev.TIMEXi,1)]
					TIMEalb <- c(min(prev.TIMEXi),max(prev.TIMEXi))
					RATEalb <- c(BASE_ALB,FINAL_ALB)
					alb.function <- approxfun(TIMEalb,RATEalb,method = "linear")
					prev.ALB <- alb.function(prev.TIMEXi)	# Linear extrapolation of ALB
					prev.ADA <- ind.bayes.data$ADA[ind.bayes.data$time == tail(prev.TIMEXi,1)]
				}
				if (covariate.scenario == "NoADA") {	# No ADA information available - assume population typical
					BASE_ALB <- ind.bayes.data$ALB[ind.bayes.data$time == prev.TIMEXi[1]]
					FINAL_ALB <- ind.bayes.data$ALB[ind.bayes.data$time == tail(prev.TIMEXi,1)]
					TIMEalb <- c(min(prev.TIMEXi),max(prev.TIMEXi))
					RATEalb <- c(BASE_ALB,FINAL_ALB)
					alb.function <- approxfun(TIMEalb,RATEalb,method = "linear")
					prev.ALB <- alb.function(prev.TIMEXi)	# Linear extrapolation of ALB
					prev.ADA <- 0
				}
				if (covariate.scenario == "NoALB") {	# No ALB information available - assume population typical
					prev.ALB <- 4
					prev.ADA <- ind.bayes.data$ADA[ind.bayes.data$time == tail(prev.TIMEXi,1)]
				}
				if (covariate.scenario == "NoCov") {	# No covariate information for ADA or ALB available
					prev.ALB <- 4
					prev.ADA <- 0
				}

			# Set up the new input data frame for Bayes estimation
				input.bayes.est.data <- data.frame(
					ID = ID.number,
					SIM = SIM.number,
					time = prev.TIMEXi,	# Time points for simulation
					ALB = prev.ALB,	# Albumin
					ADA = prev.ADA,	# Anti-drug antibodies
					ETA1 = prev.ETA1,
					ETA2 = prev.ETA2,
					ETA3 = prev.ETA3,
					ETA4 = prev.ETA4,
					ERRPRO = 0,
					amt = prev.dose,	# bayes dose
					evid = 1,	# Dosing event
					cmt = 1,	# Dose into the central compartment (compartment = 1)
					rate = -2	# Infusion duration is specific in the model file
				)
				# Make the amt given in the last time-point == 0
				# Change evid and rate accordingly
					input.bayes.est.data$amt[input.bayes.est.data$time == max(prev.TIMEXi)] <- 0
					input.bayes.est.data$evid[input.bayes.est.data$time == max(prev.TIMEXi)] <- 0
					input.bayes.est.data$rate[input.bayes.est.data$time == max(prev.TIMEXi)] <- 0
					if (max(TIMEXi) == 490) {
						input.bayes.est.data$amt[input.bayes.est.data$time == max(prev.TIMEXi)] <- prev.dose
						input.bayes.est.data$evid[input.bayes.est.data$time == max(prev.TIMEXi)] <- 1
						input.bayes.est.data$rate[input.bayes.est.data$time == max(prev.TIMEXi)] <- -2
					}
			# Modify model code ready for simulation
				prev.initial.compartment <- list(CENT = prev.init.cent,PERI = prev.init.peri,AUT = prev.init.aut)
				mod <- mod %>% init(prev.initial.compartment) %>% carry.out(SIM,amt,ERRPRO)

			# Bayesian estimation
				bayes.estimate <- function(par) {
					# Describe parameters to be optimised
						ETA1fit <- log(par[1])	# In the exponential domain
						ETA2fit <- log(par[2])
						ETA3fit <- log(par[3])
						ETA4fit <- log(par[4])
					# Add fitted parameters to the input data frame
						input.bayes.est.data$ETA1 <- ETA1fit
						input.bayes.est.data$ETA2 <- ETA2fit
						input.bayes.est.data$ETA3 <- ETA3fit
						input.bayes.est.data$ETA4 <- ETA4fit
					# Simulate concentration-time profiles with fitted parameters
						new.bayes.data <- mod %>% data_set(input.bayes.est.data) %>% mrgsim()
						new.bayes.data <- as.data.frame(new.bayes.data)
					# Pull out the predicted trough concentrations with the fitted doses for the interval
						yhat <- new.bayes.data$IPRE[new.bayes.data$time == sample.time]
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
			# Initial parameter estimates
				initial.bayes.par <- exp(c(prev.ETA1,prev.ETA2,prev.ETA3,prev.ETA4))
				par <- initial.bayes.par
				bayes.result <- optim(par,bayes.estimate,hessian = FALSE,method = "L-BFGS-B",lower = c(0.001,0.001,0.001,0.001),upper = c(Inf,Inf,Inf,Inf),control = list(parscale = par,factr = 1e7))

			# Input estimated individual parameter values and covariate values for simulation of previous interval
				input.bayes.data <- input.bayes.est.data
				new.ETA1 <- log(bayes.result$par[1])
				new.ETA2 <- log(bayes.result$par[2])
				new.ETA3 <- log(bayes.result$par[3])
				new.ETA4 <- log(bayes.result$par[4])
				input.bayes.data$ETA1 <- new.ETA1	# ETA for clearance
				input.bayes.data$ETA2 <- new.ETA2	# ETA for V1
				input.bayes.data$ETA3 <- new.ETA3	# ETA for Q
				input.bayes.data$ETA4 <- new.ETA4	# ETA for V2
			# Simulate previous interval according to the individual parameter estimates and level of covariate information
				new.bayes.data <- mod %>% data_set(input.bayes.data) %>% mrgsim()
				new.bayes.data <- as.data.frame(new.bayes.data)

			##############
			##_OPTIMISE_##
			##############
			# Optimise dose for the individual using:
				# 1. Individual BAYES predicted concentration (compartment amounts)
				# 2. Carried forward covariate information from the end of the previous interval

			# Set up a data frame to be input for dose optimisation
				input.optimise.data <- data.frame(
					ID = ID.number,
					SIM = SIM.number,
					time = TIMEX,	# Time-points for Simulation
					ALB = tail(prev.ALB,1),	# Carry forward last albumin concentration
					ADA = prev.ADA,	# Carry forward last ADA status
					ETA1 = new.ETA1,
					ETA2 = new.ETA2,
					ETA3 = new.ETA3,
					ETA4 = new.ETA4,
					ERRPRO = 0,
					amt = prev.dose,	# Dose to be optimised
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

			# Pull the amount in the compartments at the end of the previous interval
				prev.cent <- new.bayes.data$CENT[new.bayes.data$time == sample.time]
				prev.peri <- new.bayes.data$PERI[new.bayes.data$time == sample.time]
				prev.aut <- new.bayes.data$AUT[new.bayes.data$time == sample.time]
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
					# Calculate the new dose for the next interval based on "sample" and "dose"
						if (prev.DV < trough.target | prev.DV >= trough.upper) {
							initial.dose <- trough.target/prev.DV*prev.dose	# Adjust the dose if out of range
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

			# Create a data frame ready for simulating what happened to the patient given the Bayesian guided doses
				# Input optimised doses ready for simulation
					input.optimise.data$amt[input.optimise.data$time == TIMEXi[1]] <- optimised.doses$par[1]	# First dose
					input.optimise.data$amt[input.optimise.data$time == TIMEXi[2]] <- optimised.doses$par[2]	# Second dose
					input.optimise.data$amt[input.optimise.data$time == TIMEXi[3]] <- optimised.doses$par[3]	# Third dose
				# Input changing ETA values and covariate values
					ind.data <- pop.data.import[pop.data.import$ID == ID.number & pop.data.import$SIM == SIM.number & pop.data.import$TIME %in% TIMEX,]
					input.optimise.data$ETA1 <- ind.data$ETA1
					input.optimise.data$ETA2 <- ind.data$ETA2
					input.optimise.data$ETA3 <- ind.data$ETA3
					input.optimise.data$ETA4 <- ind.data$ETA4
					input.optimise.data$ALB <- ind.data$ALB
					input.optimise.data$ADA <- ind.data$ADA
					input.optimise.data$ERRPRO <- ind.data$ERRPRO
			# Simulate
				new.optimise.data <- mod %>% data_set(input.optimise.data) %>% mrgsim()
				new.optimise.data <- as.data.frame(new.optimise.data)
		}
		# Simulate concentration-time profiles for individuals in ID.data
			new.optimise.bayes.data <- ddply(ID.data, .(SIM,ID), population.bayes, .parallel = TRUE, .progress = "text")
			new.optimise.bayes.data
	}

# Simulate intervals separately
# Simulate the second interval
	optimise.bayes.data2 <- interval.bayes(TIME1i,TIME2,TIME2i,sample.time1)
# Simulate the third interval
	optimise.bayes.data3 <- interval.bayes(TIME2i,TIME3,TIME3i,sample.time2)
# Simulate the fourth interval
	optimise.bayes.data4 <- interval.bayes(TIME3i,TIME4,TIME4i,sample.time3)

# Combine bayes.dataX
	optimise.bayes.data <- rbind(conc.data,optimise.bayes.data2,optimise.bayes.data3,optimise.bayes.data4)
	optimise.bayes.data <- optimise.bayes.data[with(optimise.bayes.data, order(optimise.bayes.data$ID,optimise.bayes.data$SIM)), ]	# Sort by ID then SIM

# # ------------------------------------------------------------------------------
# # Test plot
# 	plotobj5 <- NULL
# 	plotobj5 <- ggplot(optimise.bayes.data)
# 	plotobj5 <- plotobj5 + stat_summary(aes(x = time,y = IPRE),geom = "line",fun.y = median,colour = "red")
# 	plotobj5 <- plotobj5 + stat_summary(aes(x = time,y = IPRE),geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "red",alpha = 0.3)
# 	plotobj5 <- plotobj5 + geom_hline(aes(yintercept = trough.target),linetype = "dashed")
# 	plotobj5 <- plotobj5 + geom_hline(aes(yintercept = trough.upper),linetype = "dashed")
# 	plotobj5 <- plotobj5 + scale_y_log10("Infliximab Concentration (mg/L)\n",breaks = c(0.001,0.01,0.1,1,10,100,100),labels = c(0.001,0.01,0.1,1,10,100,100))
# 	plotobj5 <- plotobj5 + scale_x_continuous("\nTime (days)")
# 	plotobj5
