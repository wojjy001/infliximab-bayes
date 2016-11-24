# Simulating Population
# R Script for simulating infliximab concentrations for a population
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list=ls(all=TRUE))

# Load package libraries
	library(ggplot2)	# Plotting
	library(grid)	# Plotting

# ------------------------------------------------------------------------------
# Set a directory for where plots can be saved (best where this R script is saved)
	setwd("/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/")

# Define a custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 16))

# Function for calculating the median, and 2.5th and 97.5th percentiles for plotting simulation results
	sumfuncx <- function(x) {
		stat1 <- median(x)
		stat2 <- quantile(x,probs = 0.025,names = F)
		stat3 <- quantile(x,probs = 0.975,names = F)
		stat4 <- length(x)
		result <- c("median" = stat1,"low" = stat2,"hi" = stat3,"n" = stat4)
		result
	}

# Set seed for reproducible results
	set.seed(123456)

# ------------------------------------------------------------------------------
# Create a time function
	TIME <- seq(from = 0,to = 546,by = 1)
# Set beginning and final albumin levels
	BASE_ALB <- 4
	FINAL_ALB <- 3
# Set specifications for first sine wave
	AMP_ALB1 <- 0.05
	FREQ_ALB1 <- 1/100
	PHASE_ALB1 <- 0
# Calculate albumin over time using the first sine wave
	TIMEalb <- c(0,546)
	RATEalb <- c(BASE_ALB,FINAL_ALB)
	step.alb <- approxfun(TIMEalb,RATEalb,method = "linear")
	ALB1 <- step.alb(TIME)*(1+AMP_ALB1*sin(2*pi*FREQ_ALB1*TIME+PHASE_ALB1))
# Plot first sine wave
	plotobj1 <- NULL
	plotobj1 <- ggplot()
	plotobj1 <- plotobj1 + geom_line(aes(x = TIME,y = ALB1),colour = "red")
	plotobj1 <- plotobj1 + scale_y_continuous(lim = c(0,6))
	print(plotobj1)

# Set specifications for second sine wave
	AMP_ALB2 <- 0.05
	FREQ_ALB2 <- 1/60
	PHASE_ALB2 <- 7
	ALB2 <- step.alb(TIME)*(1+AMP_ALB2*sin(2*pi*FREQ_ALB2*TIME+PHASE_ALB2))

	plotobj2 <- plotobj1
	plotobj2 <- plotobj2 + geom_line(aes(x = TIME,y = ALB2),colour = "blue")
	print(plotobj2)

# Add a third sine wave
	AMP_ALB3 <- 0.05
	FREQ_ALB3 <- 1/5
	PHASE_ALB3 <- 0
	ALB3 <- step.alb(TIME)*(1+AMP_ALB3*sin(2*pi*FREQ_ALB3*TIME+PHASE_ALB3))

	plot.breaks <- c(0,14,42,98,154,210,266,322,378,434,490,546)

	plotobj3 <- plotobj2
	plotobj3 <- plotobj3 + geom_line(aes(x = TIME,y = ALB3),colour = "darkgreen")
	plotobj3 <- plotobj3 + scale_y_continuous(lim = c(0,6))
	print(plotobj3)

# Add them all together
	ALB4 <- step.alb(TIME)*(1+AMP_ALB1*sin(2*pi*FREQ_ALB1*TIME+PHASE_ALB1)+AMP_ALB2*sin(2*pi*FREQ_ALB2*TIME+PHASE_ALB2)+AMP_ALB3*sin(2*pi*FREQ_ALB3*TIME+PHASE_ALB3))

	plotobj4 <- NULL
	plotobj4 <- ggplot()
	plotobj4 <- plotobj4 + geom_line(aes(x = TIME,y = ALB4))
	plotobj4 <- plotobj4 + scale_y_continuous(lim = c(0,6))
	print(plotobj4)
