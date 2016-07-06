# Time-weighted Bayes project
# Script for reading saved output data and summarising/plotting results
# ------------------------------------------------------------------------------
# Load package libraries
	library(ggplot2)	# Plotting package
	library(grid)	# Plotting package
	library(plyr)	# Split and rearrange data, ddply functions
# Custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 14))

# ------------------------------------------------------------------------------
# Set working directory
	n <- 12	# Number of seed individuals
	nsim <- 10	# Number of simulations of each individuals

	work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/"
	project.dir <- paste0(work.dir,"SIM",nsim,"_IND",n)
	setwd(project.dir)

# ------------------------------------------------------------------------------
# Functions for calculating 95% prediction intervals
	CI95lo <- function(x) quantile(x,probs = 0.025)
	CI95hi <- function(x) quantile(x,probs = 0.975)

# Summary function for calculating median and 95% prediction intervals
	summary.function95 <- function(x) {
		median <- median(x)	# median
		stat.lo95 <- CI95lo(x)	# 2.5th percentile
		stat.hi95 <- CI95hi(x)	# 97.5th percentile
		result <- c(median,stat.lo95,stat.hi95)
		names(result)[c(1,2,3)] <- c("median","stat.lo95","stat.hi95")
		result
	}

# Last per ID function
	lastperID <- function(x) tail(x,1)

# ------------------------------------------------------------------------------
###########
##_LABEL_##
###########
# Read in label_simulation.csv
	label.data <- read.csv("label_simulation.csv")

# Plot population predicted and 95% prediction intervals
# Facet for each "individual"
	plotobj1 <- NULL
	plotobj1 <- ggplot(label.data)
	plotobj1 <- plotobj1 + geom_line(aes(x = time,y = IPRE),data = label.data[label.data$SIM == 0,],colour = "red")
	plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPRE),data = label.data[label.data$SIM != 0,],geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "red",alpha = 0.3)
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 3),linetype = "dashed")
	plotobj1 <- plotobj1 + scale_x_continuous("\nTime (days)")
	plotobj1 <- plotobj1 + scale_y_log10("Infliximab Concentration (mg/L)\n")
	plotobj1 <- plotobj1 + facet_wrap(~ID)
	print(plotobj1)

# Summarise AUT over time
	label.AUT.summary <- ddply(label.data[label.data$SIM != 0,], .(time), function(label.data) summary.function95(label.data$AUT))

# Plot AUT summary data over time
	plotobj2 <- NULL
	plotobj2 <- ggplot(label.AUT.summary)
	plotobj2 <- plotobj2 + geom_line(aes(x = time,y = median),colour = "red")
	plotobj2 <- plotobj2 + geom_ribbon(aes(x = time,ymin = stat.lo95,ymax = stat.hi95),fill = "red",alpha = 0.3)
	plotobj2 <- plotobj2 + scale_x_continuous("\nTime since first dose (days)")
	plotobj2 <- plotobj2 + scale_y_continuous("Cumulative time under target trough (days)\n")
	# plotobj2 <- plotobj2 + facet_wrap(~ID)
	print(plotobj2)

#	Find final AUT for each individual
	label.data.last <- ddply(label.data, .(SIM,ID), lastperID)

# Plot boxplots
	plotobj3 <- NULL
	plotobj3 <- ggplot(label.data.last)
	plotobj3 <- plotobj3 + geom_boxplot(aes(x = ID,y = AUT,group = ID))
	plotobj3 <- plotobj3 + scale_x_continuous(breaks = unique(label.data.last$ID))
	plotobj3
