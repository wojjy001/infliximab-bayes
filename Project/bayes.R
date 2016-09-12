# in silico infliximab dosing project
# Script for estimating individual parameter values and using them to guide dose optimisation
# Subsequent doses are dependent on measured trough concentrations
# Doses are optimised using maximum likelihood estimation
# ------------------------------------------------------------------------------
# Optimise doses using maximum likelihood estimation
	# first.int.data <- first.int.data[first.int.data$SIM == 70 & first.int.data$ID == 1,]
	bayes.function <- function(first.int.data) {
		# Set up a loop that will sample the individual's concentration, estimate empirical Bayes parameters, optimise their dose and administer until time = 546 days
			# Make all predicted concentrations (IPRE) and PK parameter values after sample.time1 == NA
				conc.data <- as.data.frame(first.int.data)
				conc.data$IPRE[conc.data$time > max(sample.times)] <- NA
		# Define a variable telling the loop if previous bayes results are present
			# If they are present, they will be used as initial estimates for the next dosing interval
				previous.bayes.results.present <- FALSE
				bayes.result.list <- list("Initial" = NA)

		# If the last predicted concentration in the data frame (i.e., when time = 546) is NA, then continue with the loop
			repeat {

				###########
				##_BAYES_##
				###########
				# Estimate individual parameter values using:
					# 1. Trough sample from the end of the previous interval
					# 2. Covariate values measured at beginning and end of previous interval
					# 3. Known doses that were administered during the previous interval

					# Time of most recent sample
						last.sample <- max(sample.times)
					# Previous time-interval
						prev.TIME <- seq(from = 0,to = last.sample,by = time.int)
					# Pull covariate information from the beginning and end of the previous interval
						# Weight
							WT1 <- conc.data$WTCOV[conc.data$time %in% sample.times]
							TIMEwt <- c(sample.times)
							RATEwt <- c(WT1)
							extrap.wt <- approxfun(TIMEwt,RATEwt,method = "linear")
							prev.WT <- extrap.wt(prev.TIME)
						# Albumin
							ALB1 <- conc.data$ALBCOV[conc.data$time %in% sample.times]
							TIMEalb <- c(sample.times)
							RATEalb <- c(ALB1)
							extrap.alb <- approxfun(TIMEalb,RATEalb,method = "linear")
							prev.ALB <- extrap.alb(prev.TIME)
						# ADA status
							ADA1 <- conc.data$ADA[conc.data$time %in% sample.times]
							ada.data <- data.frame(sample.times,ADA1)
							if (tail(ADA1,1) == 1) {
								first.ada.time <- head(ada.data$sample.times[ada.data$ADA1 == 1],1)
								carry.back.to.time <- tail(ada.data$sample.times[ada.data$ADA1 == 0],1)
								ada.data$ADA1[ada.data$sample.times == carry.back.to.time] <- 1
							}
							TIMEada <- c(ada.data$sample.times)
							RATEada <- c(ada.data$ADA1)
							extrap.ada <- approxfun(TIMEada,RATEada,method = "constant")
							prev.ADA <- extrap.ada(prev.TIME)
					# Set up a new input data frame for Bayes estimation
						input.bayes.data <- conc.data[conc.data$time %in% prev.TIME,]
						input.bayes.data$TIME_WT <- prev.WT
						input.bayes.data$TIME_ADA <- prev.ADA
						input.bayes.data$TIME_ALB <- prev.ALB
						# Reduce data frame to only include sample times and dosing times
							# input.bayes.data <- input.bayes.data[input.bayes.data$time %in% sample.times | input.bayes.data$amt != 0,]
							# Re-add evid and rate columns
								input.bayes.data$cmt <- 1	# Signifies which compartment the dose goes into
								input.bayes.data$evid <- 1	# Signifies dosing event
								input.bayes.data$evid[input.bayes.data$amt == 0] <- 0
								input.bayes.data$rate <- -2	# Signifies that infusion duration is specified in model file
								input.bayes.data$rate[input.bayes.data$amt == 0] <- 0
					# Flag that this is bayes scenario
						input.bayes.data$FLAG <- 1

					# repeat {
						# Initial estimates for Bayes parameters
							initial.bayes.par <- exp(c(0,0,0,0))	# Population typical the first time bayes parameters are estimated
							# if (previous.bayes.results.present == TRUE) initial.bayes.par <- c(bayes.result$par[1],bayes.result$par[2],bayes.result$par[3],bayes.result$par[4])*exp(runif(4,min = -0.1,max = 0.1))	# Use previous parameter results as initial estimates mulitpled by a random number so that they are not exactly the same
							par <- initial.bayes.par
						# Previous DV
							prev.DV <- input.bayes.data$DV[input.bayes.data$time %in% sample.times[sample.times != 0]]
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
									new.bayes.data <- mod %>% mrgsim(data = input.bayes.data) %>% as.tbl
									new.bayes.data$IPRE[is.finite(new.bayes.data$IPRE) == F | new.bayes.data$IPRE < 1e-8] <- 1e-8
								# Pull out the predicted trough concentrations with the fitted doses for the interval
									yhat <- new.bayes.data$IPRE[new.bayes.data$time %in% sample.times[sample.times != 0]]
									# Posterior log-likelihood
									# Error model: Y = IPRE*(1+ERRPRO), Y = IPRE + IPRE*ERRPRO
										loglikpost.sd <- ERRPRO	# No time-weighting
										loglikpost <- dnorm(prev.DV,mean = yhat,sd = yhat*loglikpost.sd,log = T)
									# Prior log-likelihood
										ETA <- c(ETA1fit,ETA2fit,ETA3fit,ETA4fit)
										ETABSV <- c(PPVCL,PPVV1,PPVQ,PPVV2)
										loglikprior <- dnorm(ETA,mean = 0,sd = ETABSV,log = T)
								# Objective function value and minimise the value
									objective <- -1*sum(loglikpost,loglikprior)
							}
						# Run bayes.estimate function through optim
							bayes.result <- optim(par,bayes.estimate,hessian = FALSE,
								method = "L-BFGS-B",
								lower = c(0.001,0.001,0.001,0.001),upper = c(Inf,Inf,Inf,Inf),
								control = list(parscale = par,fnscale = bayes.estimate(par),factr = 1e12)#,
								# gr = gradient.function
							)
						# Check if bayes.result exists - may not if there were errors during optimisation
							previous.bayes.results.present <- TRUE
							bayes.result.list <- list(bayes.result.list,bayes.result)
							# Calculate concentrations according to new Bayes estimates
								# Convert new ETA values (estimated in the exp() domain)
									new.ETA1 <- log(bayes.result$par[1])
									new.ETA2 <- log(bayes.result$par[2])
									new.ETA3 <- log(bayes.result$par[3])
									new.ETA4 <- log(bayes.result$par[4])
							# Add to "input.bayes.data" data frame
								input.bayes.data$ETA1 <- new.ETA1
								input.bayes.data$ETA2 <- new.ETA2
								input.bayes.data$ETA3 <- new.ETA3
								input.bayes.data$ETA4 <- new.ETA4
							# Simulate previous interval according to the individual parameter estimates and level of covariate information
								bayes.sim.data <- mod %>% mrgsim(data = input.bayes.data,carry.out = c("amt","ERRPRO")) %>% as.tbl

				##############
				##_OPTIMISE_##
				##############
				# Optimise dose for the individual using:
					# 1. Individual BAYES predicted concentration (compartment amounts)
					# 2. Carried forward covariate information from the end of the previous interval

				# Create an input data frame for simulation
					input.sim.data <- conc.data
				# Only optimise doses if the last trough concentration from bayes.sim.data is out of target range
					if (bayes.sim.data$IPRE[bayes.sim.data$time == last.sample] < trough.target | bayes.sim.data$IPRE[bayes.sim.data$time == last.sample] >= trough.upper) {
						#	Next time-interval (default 56 days after the last sample - but dose.int can be updated via next.dose.int)
						# The next trough is the last sample time plus dosing interval time
							next.trough <- last.sample + dose.int
							next.TIME <- c(last.sample,next.trough)
							# next.TIME[next.TIME > 546] <- 546
							next.TIME.int <- seq(from = next.TIME[1],to = next.TIME[2],by = time.int)
						# Set up a data frame to be input for dose optimisation based on Bayes parameters
							input.optimise.data <- conc.data
							input.optimise.data$TIME_WT <- tail(prev.WT,1)
							input.optimise.data$TIME_ALB <- tail(prev.ALB,1)	# Carried forward last value of the previous interval
							input.optimise.data$TIME_ADA <- tail(prev.ADA,1)
							input.optimise.data$ETA1 <- new.ETA1
							input.optimise.data$ETA2 <- new.ETA2
							input.optimise.data$ETA3 <- new.ETA3
							input.optimise.data$ETA4 <- new.ETA4
						# Modify model code ready for simulation
							prev.bayes.CENT <- bayes.sim.data$CENT[bayes.sim.data$time == last.sample]	# Last value for bayes predicted CENT
							prev.bayes.PERI <- bayes.sim.data$PERI[bayes.sim.data$time == last.sample]
							prev.bayes.TUT <- bayes.sim.data$TUT[bayes.sim.data$time == last.sample]
							prev.bayes.AUT <- bayes.sim.data$AUT[bayes.sim.data$time == last.sample]
							initial.bayes.compartment <- list(CENT = prev.bayes.CENT,PERI = prev.bayes.PERI,TUT = prev.bayes.TUT,AUT = prev.bayes.AUT)
							optim.mod1 <- mod %>% init(initial.bayes.compartment) %>% carry.out(amt,ERRPRO)
						# Initial dose and error estimates
							initial.dose <- amt.min*tail(prev.WT,1)	# 5 mg/kg - label recommendation
							initial.error <- 0.01
							par <- c(initial.dose,initial.error)
						# Re-add evid and rate columns
							input.optimise.data$amt <- 0
							input.optimise.data$amt[input.optimise.data$time == last.sample] <- initial.dose
							input.optimise.data$cmt <- 1	# Signifies which compartment the dose goes into
							input.optimise.data$evid <- 1	# Signifies dosing event
							input.optimise.data$evid[input.optimise.data$amt == 0] <- 0
							input.optimise.data$rate <- -2	# Signifies that infusion duration is specified in model file
							input.optimise.data$rate[input.optimise.data$amt == 0] <- 0
						# Flag that this is bayes scenario
							input.optimise.data$FLAG <- 1
						# Subset times for only the next interval for dose optimisation
							input.optimise.data <- input.optimise.data[input.optimise.data$time %in% next.TIME,]
						# Find the doses that maximum the likelihood of trough concentrations being the target
							optimise.dose <- function(par) {
								# Add fitted parameters to the input data frame
									input.optimise.data$amt[input.optimise.data$evid == 1] <- par[1]
									err <- par[2]
								# Simulate concentration-time profiles with fitted doses
									new.optimise.data <- optim.mod1 %>% mrgsim(data = input.optimise.data) %>% as.tbl
									if (exists("new.optimise.data") == FALSE) break
									new.optimise.data$IPRE[is.finite(new.optimise.data$IPRE) == F | new.optimise.data$IPRE < 1e-8] <- 1e-8
								# Pull out the predicted trough concentrations with the fitted doses for the interval
									yhat <- new.optimise.data$IPRE[new.optimise.data$time == max(next.TIME)]
									res <- dnorm(trough.target,yhat,yhat*err,log = T)	# Minimise the error between target trough (3 mg/L) and predicted trough concentrations
								# Objective function value and minimise the value
									objective <- -1*sum(res)
							}
							optimised.doses <- optim(par,optimise.dose,hessian = FALSE,method = "L-BFGS-B",lower = c(amt.min*tail(prev.WT,1),0.0001),upper = c(amt.max*tail(prev.WT,1),Inf),control = list(parscale = par,factr = 1e12))
						# If optimised dose is reaching 50 mg/kg (max) then individual is not going to stay above target trough
						# Predict when they will hit target with the 50 mg/kg dose and then make that time the time of next dose
						# Based on Bayes parameters and carried forward covariate values!!!
							input.optimise.data$amt[input.optimise.data$time == last.sample] <- optimised.doses$par[1]
							# Simulate concentration-time profiles with fitted doses
								new.optimise.data <- optim.mod1 %>% mrgsim(data = input.optimise.data) %>% as.tbl
								new.optimise.data$IPRE[is.finite(new.optimise.data$IPRE) == F | new.optimise.data$IPRE < 1e-8] <- 1e-8
								# However, if it is predicted the individual will achieve the target trough earlier, calculate when this will be achieved using TBT (time spent under target trough)
								# Round to the nearest 7 days
									if (new.optimise.data$IPRE[new.optimise.data$time == next.trough] < trough.target) {
										TUTdiff <- diff(new.optimise.data$TUT)/7
										if (TUTdiff < 0.5) Ttarget <- 0
										if (TUTdiff >= 0.5) Ttarget <- ceiling(TUTdiff)*7	# Time under target between dose and sample
									} else {
										Ttarget <- 0
									}

								next.sample <- max(new.optimise.data$time) - Ttarget	# Time that the next concentration will be sampled at
								next.dose.int <- next.sample - new.optimise.data$time[1]	# New dosing interval for the individual based on the difference between the past and next sample times
								if (next.dose.int == 0) {	# Can't have a dosing interval of 0, therefore make the new dosing interval still every 7 days (and the dose will be optimised instead for this frequency)
									next.dose.int <- 7
									next.sample <- max(new.optimise.data$time)	# Make the time of the next concentration to be sampled as what it was originally planned
								}
						# Administer the individual the optimised dose
							input.sim.data$amt[input.sim.data$time == last.sample] <- optimised.doses$par[1]
					} else {
						# Previous dose
							prev.dose.time <- head(tail(sample.times,2),1)
							prev.dose <- conc.data$amt[conc.data$time == prev.dose.time]
							input.sim.data$amt[input.sim.data$time == last.sample] <- prev.dose
							next.sample <- last.sample + next.dose.int
					}

					# Re-add evid and rate columns
						input.sim.data$cmt <- 1	# Signifies which compartment the dose goes into
						input.sim.data$evid <- 1	# Signifies dosing event
						input.sim.data$evid[input.sim.data$amt == 0] <- 0
						input.sim.data$rate <- -2	# Signifies that infusion duration is specified in model file
						input.sim.data$rate[input.sim.data$amt == 0] <- 0
					# Flag that this is simulation and want covariates to change depending on concentrations
						input.sim.data$FLAG <- 0
					# Simulate
						conc.data <- mod %>% mrgsim(data = input.sim.data,carry.out = c("amt","ERRPRO")) %>% as.tbl
						conc.data$IPRE[is.finite(conc.data$IPRE) == F | conc.data$IPRE < 1e-8] <- 1e-8
						conc.data$DV[is.finite(conc.data$DV) == F | conc.data$DV < 1e-8] <- 1e-8
					# Add the "next.sample" time to the list of sample.times
						sample.times <- sort(c(unique(c(sample.times,next.sample))))
					# Make all predicted concentrations (IPRE) and PK parameter values after sample.time1 == NA
						conc.data$IPRE[conc.data$time > max(sample.times)] <- NA

				# Stop the loop if IPRE at the last time-point has been calculated
					if (is.na(conc.data$IPRE[conc.data$time == last.time]) == FALSE) break
			}	# Brackets closing "repeat"
		conc.data
	}	# Brackets closing "bayes.function"

	optimise.bayes.data <- ddply(first.int.data, .(SIM,ID), bayes.function, .parallel = FALSE, .progress = "text", .inform = TRUE)
	# print(bayes.result.list)
# # Numerical summary of individual
# # Infusion times and amounts
# 	print(
# 		data.frame(
# 			ID = 1,
# 			times = conc.data$time[conc.data$amt != 0 | conc.data$time == 546],
# 			amt = conc.data$amt[conc.data$amt != 0 | conc.data$time == 546],
# 			int = c(0,diff(conc.data$time[conc.data$amt != 0 | conc.data$time == 546])),
# 			ADA = conc.data$ADA[conc.data$amt != 0 | conc.data$time == 546],
# 			ALB = conc.data$ALBCOV[conc.data$amt != 0 | conc.data$time == 546],
# 			WT = conc.data$WTCOV[conc.data$amt != 0 | conc.data$time == 546]
# 		)
# 	)

# ------------------------------------------------------------------------------
# Write clinical.data to a .csv file
	optimise.bayes.data.filename <- "optimise_bayes_data.csv"
	write.csv(optimise.bayes.data,file = optimise.bayes.data.filename,na = ".",quote = F,row.names = F)
