# Simulate test population created from "create_test_population.R"
# Read in the simulated file created by NONMEM and compare to mrgsolve
# Compare the results of the two simulation methods
# Prepare simulation output ready for Bayesian estimation of individual PK parameters
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))

# Set working directory
 	# work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/NONMEM_Check/"	# Mac directory
	work.dir <- "E:/Wojciechowski/infliximab-bayes/NONMEM_Check/"	# Server directory
	setwd(work.dir)
	
# Load necessary package libraries
	library(ggplot2)
	library(grid)
# Custom ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))	

# Source and compile mrgsolve model file
	source("mrgsolve_infliximab.R")

# Set seed for reproducible results
 	set.seed(123456)
	
# Model specific information required for Bayesian estimation	
	sample.times <- c(98,322,546)	# Assign sample times
	PPVCL <- sqrt(0.106929)	# Model's PPV for clearance
	PPVV1 <- sqrt(0.0225)	# Model's PPV for volume of central compartment
	PPVQ <- sqrt(1.21)	# Model's PPV for inter-compartmental clearance
	PPVV2 <- sqrt(0.638401)	# Model's PPV for volume of peripheral compartment
	model.ERRPRO <- sqrt(0.175561)	# Model's proportion RUV		
	
# ------------------------------------------------------------------------------
# Read in .csv file of test population simulation output
	input.bayes.data <- read.csv(file = "mrgsolve_simulation_output.csv")
	input.bayes.data$cmt <- 1	# Administer doses/observe compartment 1
	input.bayes.data$evid <- 0	# Signify observation events
	input.bayes.data$evid[input.bayes.data$amt != 0] <- 1 # Signify dosing events
	input.bayes.data$rate <- 0
	input.bayes.data$rate[input.bayes.data$evid == 1] <- -2	# Rate is specified as duration in model file

# Function for Bayesian estimating individual PK parameters for each individual in the dataset
	bayes.function <- function(input.bayes.data) {		
		run.once <- FALSE
		repeat {
			# Initial estimates for Bayes parameters
				initial.bayes.par <- exp(c(0,0,0,0))	# Population typical the first time bayes parameters are estimated
				if (run.once == TRUE) initial.bayes.par <- initial.bayes.par*exp(runif(4,min = -0.01,max = 0.01))	# Use previous parameter results as initial estimates mulitpled by a random number so that they are not exactly the same
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
						new.bayes.data$IPRE[is.finite(new.bayes.data$IPRE) == F | new.bayes.data$IPRE < .Machine$double.eps] <- .Machine$double.eps
					# Pull out the predicted trough concentrations with the fitted doses for the interval
						yhat <- new.bayes.data$IPRE[new.bayes.data$time %in% sample.times[sample.times != 0]]
						# Posterior log-likelihood
						# Error model: Y = IPRE*(1+ERRPRO), Y = IPRE + IPRE*ERRPRO
							loglikpost.sd <- model.ERRPRO	# No time-weighting
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
					control = list(parscale = par,fnscale = bayes.estimate(par))#,
					# gr = gradient.function
				)
				run.once <- TRUE
			if (bayes.result$message == "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH") break
		}	# Brackets closing "repeat"
		# Incorporate results into a data frame for ddply to bind
			result <- data.frame(
				ETA_CL = log(bayes.result$par[1]),
				ETA_V1 = log(bayes.result$par[2]),
				ETA_Q = log(bayes.result$par[3]),
				ETA_V2 = log(bayes.result$par[4]),
				OFV = bayes.result$value[1],
				message = bayes.result$message[1]
			)
	}
	
# Apply the bayes.function to each individual	
	r.bayes.data <- ddply(input.bayes.data, .(ID), bayes.function, .progress = "text")
	write.csv(r.bayes.data,file = paste0(work.dir,"mrgsolve_bayes_output.csv"),na = ".",quote = F,row.names = F) # Write to .csv file

# ------------------------------------------------------------------------------
# Read in the processed estimated data from NONMEM
	nonmem.dir <- paste0(work.dir,"nonmem_bayes_infliximab.nm7/")
	nonmem.bayes.data <- read.csv(paste0(nonmem.dir,"nonmem_bayes_infliximab.fit.csv"))
	nonmem.bayes.data <- nonmem.bayes.data[nonmem.bayes.data$TIME %in% sample.times & nonmem.bayes.data$AMT == 0,]
	nonmem.bayes.data <- nonmem.bayes.data[,c(1,2,11,15,17,19,21,23)] # Remove unnecessary columns
	
# Calculate the objective function value given the parameters from NONMEM
	calc.ofv.function <- function(input) {
		loglikpost <- dnorm(input$DV,mean = input$IPRE,sd = input$IPRE*model.ERRPRO,log = T)
		ETA <- c(input$ETA1[1],input$ETA2[1],input$ETA3[1],input$ETA4[1])
		ETABSV <- c(PPVCL,PPVV1,PPVQ,PPVV2)
		loglikprior <- dnorm(ETA,mean = 0,sd = ETABSV,log = T)
		objective <- -1*sum(loglikpost,loglikprior)
		result <- data.frame(
			ETA_CL = input$ETA1[1],
			ETA_V1 = input$ETA2[1],
			ETA_Q = input$ETA3[1],
			ETA_V2 = input$ETA4[1],
			OFV = objective
		)
	}
	nonmem.bayes.data <- ddply(nonmem.bayes.data, .(ID), calc.ofv.function)
	
# ------------------------------------------------------------------------------
# Compare the R and NONMEM results
# Calculate absolute error respect to NONMEM results
# Has to be absolute, because values cross over zero
	ETA_CL.err <- r.bayes.data$ETA_CL-nonmem.bayes.data$ETA_CL
	ETA_V1.err <- r.bayes.data$ETA_V1-nonmem.bayes.data$ETA_V1
	ETA_Q.err <- r.bayes.data$ETA_Q-nonmem.bayes.data$ETA_Q
	ETA_V2.err <- r.bayes.data$ETA_V2-nonmem.bayes.data$ETA_V2
	OFV.err <- r.bayes.data$OFV-nonmem.bayes.data$OFV
	ETA.err <- list("ETA_CL.err" = ETA_CL.err,"ETA_V1.err" = ETA_V1.err,"ETA_Q.err" = ETA_Q.err,"ETA_V2.err" = ETA_V2.err,"OFV.err" = OFV.err)
# Calculate the median and 95% confidence intervals
	summary.function <- function(input) {
		median.value <- median(input)
		CI50lo.value <- quantile(input,probs = 0.25,names = F)
		CI50hi.value <- quantile(input,probs = 0.75,names = F)
		summary.result <- c(median.value,CI50lo.value,CI50hi.value)
		names(summary.result) <- c("median","CI50lo","CI50hi")
		summary.result
	}
	summary.result <- ldply(ETA.err,summary.function)

# Plot the results
	plotobj <- NULL
	plotobj <- ggplot(summary.result)
	plotobj <- plotobj + geom_point(aes(x = .id,y = median,colour = .id),size = 4)
	plotobj <- plotobj + geom_errorbar(aes(x = .id,ymin = CI50lo,ymax = CI50hi,colour = .id),size = 1)
	plotobj <- plotobj + geom_hline(aes(yintercept = 0),linetype = "dashed")
	plotobj <- plotobj + scale_y_continuous("Absolute Error of R with Respect to NONMEM\n",lim = c(-0.05,0.05))
	plotobj <- plotobj + scale_x_discrete("\nParameter")
	plotobj <- plotobj + theme(legend.position = "none")
	plotobj

  ggsave(plot = plotobj,filename = paste0(work.dir,"bayes_estimates_error.png"),units = "cm",width = 20,height = 20,dpi = 300)