# Dose optimisation script
# In my time-weighted Bayes project, there are some people with such low trough concentrations that they are bugging out the dose optimisation algorithm
# I need to find out why this is happening and fix it!
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
	input.optimise.data <- data.frame(
		ID = 1,
		time = c(378,434,490,546),	# Infusion/trough times
		ALB = 3.621317186,	# Albumin (g/dL)
		ADA = 1,	# Anti-drug antibodies
		ETA1 = 0.688112345,	# ETA for CL
		ETA2 = -0.087150995,	# ETA for V1
		ETA3 = 0.764270184,	# ETA for Q
		ETA4 = -0.964647963,	# ETA for V2
		ERRPRO = -1.03772694,	# Random effect for residual error
		amt = 0,	# Dose to be optimised
		evid = c(1,1,1,0),	# Dosing event
		cmt = 1,	# Dose into the central compartment (compartment = 1)
		rate = c(-2,-2,-2,0)	# Infusion duration is specified in the model file
	)

# Initial compartment values
	initial.compartment <- list(CENT = 0.53246085,PERI = 0.124347732,AUT = 84.23292059)
	optim.mod1 <- mod %>% init(initial.compartment) %>% carry.out(amt,ERRPRO)

# Intiial parameter estimates for dose optimisation
	last.dose <- 468135.7002
	last.IPRE <- 0.18832934
	# Calculate the new dose for the next interval based on "sample" and "dose"
		if (last.IPRE < 3 | last.IPRE >= 5) {
			initial.dose <- 3/last.IPRE*last.dose	# Adjust the dose if out of range
		} else {
			initial.dose <- last.dose	# Continue with previous dose if within range
		}
		if (initial.dose > 1000000) {
			initial.dose <- NA
		}

	if (is.na(initial.dose) == FALSE & initial.dose != last.dose) {
		# Limits of parameters
			initial.par <- c(initial.dose,initial.dose,initial.dose,0.01)
			lower.limit <- c(0.0001,0.0001,0.0001,0.0001)
			upper.limit <- c(Inf,Inf,Inf,Inf)
			par <- initial.par

		# Find the doses that maximum the likelihood of trough concentrations being the target
			optimise.dose <- function(par) {
			# Add fitted parameters to the input data frame
				input.optimise.data$amt[1] <- par[1]
				input.optimise.data$amt[2] <- par[2]
				input.optimise.data$amt[3] <- par[3]
				err <- par[4]
			# Simulate concentration-time profiles with fitted doses
				new.optimise.data <- optim.mod1 %>% data_set(input.optimise.data) %>% mrgsim()
				new.optimise.data <- as.data.frame(new.optimise.data)
			# Pull out the predicted trough concentrations with the fitted doses for the interval
				yhat <- new.optimise.data$IPRE[new.optimise.data$time > 378]
				print(yhat)
				res <- dnorm(3,yhat,yhat*err,log = T)	# Minimise the error between target trough (3 mg/L) and predicted trough concentrations
			# Objective function value and minimise the value
				objective <- -1*sum(res)

				diagnosticflag <- T

			    if (diagnosticflag == T) {
			    # Plot parameter values during fitting
			      bayes1buffer <<- c(bayes1buffer,log10(par[1]/initial.par[1]))
			      bayes2buffer <<- c(bayes2buffer,log10(par[2]/initial.par[2]))
						bayes3buffer <<- c(bayes3buffer,log10(par[3]/initial.par[3]))
						bayes4buffer <<- c(bayes4buffer,log10(par[4]/initial.par[4]))
			      iterationbuffer <- 1:length(bayes1buffer)
			      if (length(iterationbuffer) %% 5 == 0) { # Plot every 5th iteration
		          plot(iterationbuffer,bayes1buffer,type = 'l',xlim = c(1,length(iterationbuffer)),ylim = c(-5,5),col = "red",xlab = "Iteration",ylab = "Normalised Parameter Value",main = "Estimating Bayes Parameters\n")
		          points(iterationbuffer,bayes2buffer,type = 'l',col = "blue")
							points(iterationbuffer,bayes3buffer,type = 'l',col = "darkgreen")
							points(iterationbuffer,bayes4buffer,type = 'l',col = "orange")
		          legend(1,-1,legend = c("DOSE1","DOSE2","DOSE3","ERR"),lty = c(1,1),col = c("red","blue","darkgreen","orange"))
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

			optimised.doses <- optim(par,optimise.dose,hessian = FALSE,method = "L-BFGS-B",lower = lower.limit,upper = upper.limit,control = list(parscale = par,factr = 1e12))
	} else {
		optimised.doses <- NA
	}

	# Extract doses from "optimised.doses"
		if (is.na(initial.dose) == FALSE & initial.dose != last.dose) {
			DOSE1 <- optimised.doses$par[1]
			DOSE2 <- optimised.doses$par[2]
			if (interval == 2) {
				DOSE3 <- NA
				ERR <- optimised.doses$par[3]
			}
			if (interval != 2) {	# Other intervals contain three doses
				DOSE3 <- optimised.doses$par[3]
				ERR <- optimised.doses$par[4]
			}
		} else if (is.na(initial.dose) == FALSE & initial.dose == last.dose) {
			DOSE1 <- last.dose
			DOSE2 <- last.dose
			DOSE3 <- NA
			if (interval != 2) {	# Other intervals contain three doses
				DOSE3 <- last.dose
			}
			ERR <- 0
		} else {
			DOSE1 <- 1000000
			DOSE2 <- 1000000
			DOSE3 <- NA
			if (interval != 2) {	# Other intervals contain three doses
				DOSE3 <- 1000000
			}
			ERR <- 0
		}
