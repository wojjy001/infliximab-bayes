# Bayesian Estimation script
# In my time-weighted Bayes project, there are some individuals with large ETA values that are bugging out when estimating Bayes parameters
# May be an issue with initial estimates of parameters for these individuals
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list=ls(all=TRUE))

# Load package libraries
	library(plyr) # Split and rearrange data, ddply function
	library(dplyr)	# Split and rearrange data
	library(mrgsolve)	# Metrum Research Group differential equation solver

# Set working directory
	setwd("/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/")

# ------------------------------------------------------------------------------
# Source the infliximab model file
	source("/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/model.R")

# Set up an input data frame for optimisation algorithm and differential equation solver
	input.bayes.data <- data.frame(
		ID = 1,
		time = c(0,14,42,98,154,210),	# Infusion/trough times
		ALB = 4,	# Albumin (g/dL)
		ADA = 0,	# Anti-drug antibodies
		ETA1 = 0,	# ETA for CL
		ETA2 = 0,	# ETA for V1
		ETA3 = 0,	# ETA for Q
		ETA4 = 0,	# ETA for V2
		ERRPRO = 0,	# Random effect for residual error
		amt = c(350,350,350,97802.66728,97789.57663,0),	# Dose to be optimised
		evid = c(1,1,1,1,1,0),	# Dosing event
		cmt = 1,	# Dose into the central compartment (compartment = 1)
		rate = c(-2,-2,-2,-2,-2,0)	# Infusion duration is specified in the model file
	)

	# Limits of parameters
		initial.bayes.par <- exp(c(0.668791497,-0.103271868,0.726685569,-1.019384554))
		lower.limit <- c(0.0001,0.0001,0.0001,0.0001)
		upper.limit <- c(Inf,Inf,Inf,Inf)
		par <- initial.bayes.par

		prev.DV <- c(0.009834309*exp(0.01829042),4.343537334*exp(0.019383009))
	# Values for PPV (Population Parameter Variability), as SDs
		PPVCL <- 0.327
		PPVV1 <- 0.150
		PPVQ <- 1.10
		PPVV2 <- 0.799
	# Value for RUV (Residual Unexplained Variability), as SD
		ERRPRO <- 0.419

# Estimate individual PK parameters
	method.scenario <- "Peck1.01"
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
			yhat <- new.bayes.data$IPRE[new.bayes.data$time %in% c(98,210)]
			# Posterior log-likelihood
			# Error model: Y = IPRE*exp(ERRPRO), log(Y) = log(IPRE) + ERRPRO
				TIMET <- max(new.bayes.data$time) - new.bayes.data$time	# Time since last observation
				if (method.scenario == "NTimeWeight") loglikpost.sd <- ERRPRO	# No time-weighting
				if (method.scenario == "Peck1.005") loglikpost.sd <- ERRPRO*1.005^TIMET  # Peck method, Q = 1.005
				if (method.scenario == "Peck1.01")	loglikpost.sd <- ERRPRO*1.01^TIMET	# Peck method, Q = 1.01
				loglikpost <- dnorm(log(prev.DV),mean = log(yhat),sd = loglikpost.sd,log = T)
			# Prior log-likelihood
				ETA <- c(ETA1fit,ETA2fit,ETA3fit,ETA4fit)
				ETABSV <- c(PPVCL,PPVV1,PPVQ,PPVV2)
				loglikprior <- dnorm(ETA,mean = 0,sd = ETABSV,log = T)
		# Objective function value and minimise the value
			objective <- -1*sum(loglikpost,loglikprior)

			diagnosticflag <- T

	    if (diagnosticflag == T) {
	    # Plot parameter values during fitting
	      bayes1buffer <<- c(bayes1buffer,log10(par[1]/initial.bayes.par[1]))
	      bayes2buffer <<- c(bayes2buffer,log10(par[2]/initial.bayes.par[2]))
				bayes3buffer <<- c(bayes3buffer,log10(par[3]/initial.bayes.par[3]))
				bayes4buffer <<- c(bayes4buffer,log10(par[4]/initial.bayes.par[4]))
	      iterationbuffer <- 1:length(bayes1buffer)
	      if (length(iterationbuffer) %% 5 == 0) { # Plot every 5th iteration
          plot(iterationbuffer,bayes1buffer,type = 'l',xlim = c(1,length(iterationbuffer)),ylim = c(-5,5),col = "red",xlab = "Iteration",ylab = "Normalised Parameter Value",main = "Estimating Bayes Parameters\n")
          points(iterationbuffer,bayes2buffer,type = 'l',col = "blue")
					points(iterationbuffer,bayes3buffer,type = 'l',col = "darkgreen")
					points(iterationbuffer,bayes4buffer,type = 'l',col = "orange")
          legend(1,-1,legend = c("ETA1","ETA2","ETA3","ETA4"),lty = c(1,1),col = c("red","blue","darkgreen","orange"))
	    	}
			}
		objective
	}

# Initial buffer values for rolling plot
  plot.new() # Clear plot for GUI
  bayes1buffer <<- NA  # Global assignment here
  bayes2buffer <<- NA
	bayes3buffer <<- NA
	bayes4buffer <<- NA

	bayes.eta <- optim(par,bayes.estimate,hessian = FALSE,method = "L-BFGS-B",lower = lower.limit,upper = upper.limit,control = list(parscale = par,factr = 1e12))
	bayes.eta
