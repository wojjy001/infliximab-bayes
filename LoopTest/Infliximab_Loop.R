# Script that wraps an individual patient into a simulating-fitting-optimisation loop
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list=ls(all=TRUE))

# Load package libraries
	library(plyr) # Split and rearrange data, ddply function
	library(dplyr)	# Split and rearrange data
	library(mrgsolve)	# Metrum Research Group differential equation solver
	library(ggplot2)	# Plotting
	library(grid)	# Plotting
	library(numDeriv)	# Package containing function that returns numerical gradients of Bayes objective function value

# Set working directory
	work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/LoopTest/"
	setwd(work.dir)

# ------------------------------------------------------------------------------
# Source the infliximab first_interval file
	source(paste0(work.dir,"first_interval.R"))
	# Output object is called "first.int.data"
	# Sample a random individual from the simulated data
		sample.ID <- sample(c(1,6,11,16),1)
		# sample.ID <- 6
		conc.data <- first.int.data[first.int.data$SIM == 1 & first.int.data$ID == sample.ID,] # First three doses for the individual have been simulated

# ------------------------------------------------------------------------------
# Set up a loop that will sample the individual's concentration, estimate empirical Bayes parameters, optimise their dose and administer until time = 546 days
# Define the last time-point to be simulated
	last.time <- 546	# days
# After the initiation phase, the first sample will be collected at day 98
	sample.times <- c(0,98)	# days
# Initial dosing interval for the maintenance phase
	dose.int <- 56	# days
	next.dose.int <- 56	# days
# Make all predicted concentrations (IPRE) and PK parameter values after sample.time1 == NA
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
				prev.TIME <- seq(from = 0,to = last.sample,by = 1)
			# Pull covariate information from the beginning and end of the previous interval
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
				# Weight
					WT1 <- conc.data$WTCOV[conc.data$time %in% sample.times]
					TIMEwt <- c(sample.times)
					RATEwt <- c(WT1)
					extrap.wt <- approxfun(TIMEwt,RATEwt,method = "linear")
					prev.WT <- extrap.wt(prev.TIME)
			# Set up a new input data frame for Bayes estimation
				input.bayes.data <- conc.data[conc.data$time %in% prev.TIME,]
				input.bayes.data$ADA <- prev.ADA
				input.bayes.data$ALBCOV <- prev.ALB
				input.bayes.data$WTCOV <- prev.WT
				# Reduce data frame to only include sample times and dosing times
					# input.bayes.data <- input.bayes.data[input.bayes.data$time %in% sample.times | input.bayes.data$amt != 0,]
					# Re-add evid and rate columns
						input.bayes.data$cmt <- 1	# Signifies which compartment the dose goes into
						input.bayes.data$evid <- 1	# Signifies dosing event
						input.bayes.data$evid[input.bayes.data$amt == 0] <- 0
						input.bayes.data$rate <- -2	# Signifies that infusion duration is specified in model file
						input.bayes.data$rate[input.bayes.data$amt == 0] <- 0

			repeat {
				# Initial estimates for Bayes parameters
					initial.bayes.par <- exp(c(0,0,0,0))	# Population typical
					if (previous.bayes.results.present == TRUE) initial.bayes.par <- c(bayes.result$par[1],bayes.result$par[2],bayes.result$par[3],bayes.result$par[4])*exp(runif(4,min = -0.1,max = 0.1))	# Use previous parameter results as initial estimates mulitpled by a random number so that they are not exactly the same
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
							new.bayes.data <- mod.bayes %>% data_set(input.bayes.data) %>% mrgsim()
							new.bayes.data <- as.data.frame(new.bayes.data)
							new.bayes.data$IPRE[is.finite(new.bayes.data$IPRE) == F | new.bayes.data$IPRE < 0.001] <- 0.001
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
				# Gradient function for "bayes.estimate"
				# Gradient function must have the same arguments as "bayes.estimate"
					gradient.function <- function(par) {
						grad(func = bayes.estimate,x = par)
					}
				# Run bayes.estimate function through optim
					bayes.result <- optim(par,bayes.estimate,hessian = FALSE,method = "L-BFGS-B",lower = c(0.001,0.001,0.001,0.001),upper = c(Inf,Inf,Inf,Inf),control = list(parscale = par,fnscale = bayes.estimate(par),factr = 1e12),gr = gradient.function)
					# bayes.result <- optim(par,bayes.estimate,hessian = FALSE)
				if (bayes.result$convergence == 0) break
			}

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
				sim.mod1 <- mod.bayes %>% carry.out(amt,SIM,ERRPRO)
				bayes.sim.data <- sim.mod1 %>% data_set(input.bayes.data) %>% mrgsim()
				bayes.sim.data <- as.data.frame(bayes.sim.data)

		##############
		##_OPTIMISE_##
		##############
		# Optimise dose for the individual using:
			# 1. Individual BAYES predicted concentration (compartment amounts)
			# 2. Carried forward covariate information from the end of the previous interval
		# OR:
			# 1. Individual IPRE predicted concentration
			# 2. Actual covariate information

		# Create an input data frame for simulation
			input.sim.data <- conc.data
		# Only optimise doses if the last trough concentration from bayes.sim.data is out of target range
			if (bayes.sim.data$IPRE[bayes.sim.data$time == last.sample] < trough.target | bayes.sim.data$IPRE[bayes.sim.data$time == last.sample] >= trough.upper) {
				#	Next time-interval (default 56 days after the last sample - but dose.int can be updated via next.dose.int)
				# The next trough is the last sample time plus dosing interval time
					next.trough <- last.sample+dose.int
					next.TIME <- c(last.sample,next.trough)
					# next.TIME[next.TIME > 546] <- 546
					next.TIME.int <- seq(from = next.TIME[1],to = next.TIME[2],by = 1)
				# Set up a data frame to be input for dose optimisation based on Bayes parameters
					input.optimise.data <- conc.data
					input.optimise.data$ADA <- tail(prev.ADA,1)
					input.optimise.data$ETA1 <- new.ETA1
					input.optimise.data$ETA2 <- new.ETA2
					input.optimise.data$ETA3 <- new.ETA3
					input.optimise.data$ETA4 <- new.ETA4
				# Modify model code ready for simulation
					prev.bayes.CENT <- bayes.sim.data$CENT[bayes.sim.data$time == last.sample]	# Last value for bayes predicted CENT
					prev.bayes.PERI <- bayes.sim.data$PERI[bayes.sim.data$time == last.sample]
					prev.bayes.TUT <- bayes.sim.data$TUT[bayes.sim.data$time == last.sample]
					prev.bayes.AUT <- bayes.sim.data$AUT[bayes.sim.data$time == last.sample]
					prev.bayes.ALB <- tail(prev.ALB,1)	# Carried forward last value of the previous interval
					prev.bayes.WT <- tail(prev.WT,1)
					initial.bayes.compartment <- list(CENT = prev.bayes.CENT,PERI = prev.bayes.PERI,TUT = prev.bayes.TUT,AUT = prev.bayes.AUT,ALB = prev.bayes.ALB,WT = prev.bayes.WT)
					optim.mod1 <- mod.bayes %>% init(initial.bayes.compartment) %>% carry.out(amt,ERRPRO)
				# Initial dose and error estimates
					initial.dose <- 350	# 5 mg/kg - label recommendation
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
				# Subset times for only the next interval for dose optimisation
					input.optimise.data <- input.optimise.data[input.optimise.data$time %in% next.TIME,]
				# Find the doses that maximum the likelihood of trough concentrations being the target
					optimise.dose <- function(par) {
						# Add fitted parameters to the input data frame
							input.optimise.data$amt[input.optimise.data$evid == 1] <- par[1]
							err <- par[2]
						# Simulate concentration-time profiles with fitted doses
							new.optimise.data <- optim.mod1 %>% data_set(input.optimise.data) %>% mrgsim()
							new.optimise.data <- as.data.frame(new.optimise.data)
						# Pull out the predicted trough concentrations with the fitted doses for the interval
							new.optimise.data$IPRE[is.finite(new.optimise.data$IPRE) == F | new.optimise.data$IPRE < 0.001] <- 0.001
							yhat <- new.optimise.data$IPRE[new.optimise.data$time == max(next.TIME)]
							res <- dnorm(trough.target,yhat,yhat*err,log = T)	# Minimise the error between target trough (3 mg/L) and predicted trough concentrations
						# Objective function value and minimise the value
							objective <- -1*sum(res)
					}
					optimised.doses <- optim(par,optimise.dose,hessian = FALSE,method = "L-BFGS-B",lower = c(0.0001,0.0001),upper = c(3500,Inf),control = list(parscale = par,factr = 1e12))
				# If optimised dose is reaching 3500 mg (max) then individual is not going to stay above target trough
				# Predict when they will hit target with the 3500 mg dose and then make that time the time of next dose
				# Based on Bayes parameters and carried forward covariate values!!!
					input.optimise.data$amt[input.optimise.data$time == last.sample] <- optimised.doses$par[1]
					# Simulate concentration-time profiles with fitted doses
						new.optimise.data <- optim.mod1 %>% data_set(input.optimise.data) %>% mrgsim()
						new.optimise.data <- as.data.frame(new.optimise.data)
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
					next.sample <- last.sample+next.dose.int
			}

			# Re-add evid and rate columns
				input.sim.data$cmt <- 1	# Signifies which compartment the dose goes into
				input.sim.data$evid <- 1	# Signifies dosing event
				input.sim.data$evid[input.sim.data$amt == 0] <- 0
				input.sim.data$rate <- -2	# Signifies that infusion duration is specified in model file
				input.sim.data$rate[input.sim.data$amt == 0] <- 0
			# Simulate
				conc.data <- mod.sim %>% carry.out(amt,SIM,ERRPRO) %>% data_set(input.sim.data) %>% mrgsim()
				conc.data <- as.data.frame(conc.data)
			# Add the "next.sample" time to the list of sample.times
				sample.times <- sort(c(unique(c(sample.times,next.sample))))
			# Make all predicted concentrations (IPRE) and PK parameter values after sample.time1 == NA
				conc.data$IPRE[conc.data$time > max(sample.times)] <- NA
			# Plot individual's concentrations as they are being administered and optimised
				scale.log10.labels <- c(0.01,0.1,1,10,100,1000)
				plotobj2 <- NULL
				plotobj2 <- ggplot()
				plotobj2 <- plotobj2 + geom_line(aes(x = time,y = IPRE),data = bayes.sim.data[bayes.sim.data$time <= last.time,],colour = "blue",linetype = "dashed")
				plotobj2 <- plotobj2 + geom_line(aes(x = time,y = IPRE),data = conc.data[conc.data$time <= last.time,],colour = "red")
				plotobj2 <- plotobj2 + geom_point(aes(x = sample.times[!c(sample.times %in% c(0,next.sample))],y = prev.DV),size = 2)
				plotobj2 <- plotobj2 + geom_hline(aes(yintercept = trough.target),linetype = "dashed")
				plotobj2 <- plotobj2 + geom_hline(aes(yintercept = trough.upper),linetype = "dashed")
				plotobj2 <- plotobj2 + scale_y_log10("Infliximab Concentration (mg/L)\n",breaks = scale.log10.labels,labels = scale.log10.labels,lim = c(0.01,10000))
				plotobj2 <- plotobj2 + scale_x_continuous("\nTime (days)",breaks = c(0,14,42,98,154,210,266,322,378,434,490,546))
				print(plotobj2)
				# ggsave(plot = plotobj2,file = paste0("infliximab_loop_",last.sample,".png"))

		# Stop the loop if IPRE at the last time-point has been calculated
			if (is.na(conc.data$IPRE[conc.data$time == last.time]) == FALSE) break

	}	# Brackets closing "repeat"

# ------------------------------------------------------------------------------
# Numerical summary of individual
# Infusion times and amounts
	print(
		data.frame(
			ID = sample.ID,
			times = conc.data$time[conc.data$amt != 0],
			amt = conc.data$amt[conc.data$amt != 0],
			int = c(0,diff(conc.data$time[conc.data$amt != 0])),
			ADA = conc.data$ADA[conc.data$amt != 0],
			ALB = conc.data$ALBCOV[conc.data$amt != 0],
			WT = conc.data$WTCOV[conc.data$amt != 0]
		)
	)
#	Cumulative time under target trough by time = 546
 	print(conc.data$TUT[conc.data$time == last.time])
# Print cumulative time under target trough by time = 98
	# Remove the effects of lame initiation phase dosing
		print(conc.data$TUT[conc.data$time == last.time]-conc.data$TUT[conc.data$time == 98])
# Print cumulative area under target trough at time = 546
	print(conc.data$AUT[conc.data$time == last.time])
	# See the effects of lame initiation phase dosing
		print(conc.data$AUT[conc.data$time == last.time]-conc.data$AUT[conc.data$time == 98])
