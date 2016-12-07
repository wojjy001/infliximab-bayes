# Simulate test population created from "create_test_population.R"
# Read in the simulated file created by NONMEM and compare to mrgsolve
# Compare the results of the two simulation methods
# Prepare simulation output ready for Bayesian estimation of individual PK parameters
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))

# Set working direcctory
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
	
# ------------------------------------------------------------------------------
# Read in .csv file of test population
	input.data <- read.csv(file = "mrgsolve_simulation_input.csv")
	
# Simulate concentrations using mrgsolve
	mrgsolve.data <- mod %>% mrgsim(data = input.data,carry.out = c("amt,ERRPRO")) %>% as.tbl
	write.csv(mrgsolve.data,file = paste0(work.dir,"mrgsolve_simulation_output.csv"),na = ".",quote = F,row.names = F) # Write to .csv file
	
# ------------------------------------------------------------------------------
# Read in the processed simulated data from NONMEM
	nonmem.dir <- paste0(work.dir,"nonmem_simulation_infliximab.nm7/")
	nonmem.data <- read.csv(paste0(nonmem.dir,"nonmem_simulation_infliximab.fit.csv"))
	
# ------------------------------------------------------------------------------
# Plot for comparisons
	# Function to calculate numerical summaries
		summary.function <- function(input) {
			median.value <- median(input)
			CI95lo.value <- quantile(input,probs = 0.025,names = F)
			CI95hi.value <- quantile(input,probs = 0.975,names = F)
			summary.result <- c(median.value,CI95lo.value,CI95hi.value)
			names(summary.result) <- c("median","CI95lo","CI95hi")
			summary.result
		}
	# Calculate the summary for each method
		mrgsolve.summary <- ddply(mrgsolve.data, .(time), function(mrgsolve.data) summary.function(mrgsolve.data$IPRE))
		nonmem.summary <- ddply(nonmem.data, .(TIME), function(nonmem.data) summary.function(nonmem.data$IPRE))
	# Plot the results
		plotobj <- NULL
		plotobj <- ggplot()
		plotobj <- plotobj + geom_ribbon(aes(x = time,ymin = CI95lo,ymax = CI95hi),data = mrgsolve.summary,fill = "red",alpha = 0.3)
		plotobj <- plotobj + geom_ribbon(aes(x = TIME,ymin = CI95lo,ymax = CI95hi),data = nonmem.summary,fill = "blue",alpha = 0.3)
		plotobj <- plotobj + geom_line(aes(x = time,y = median),data = mrgsolve.summary,colour = "red",size = 1)
		plotobj <- plotobj + geom_line(aes(x = TIME,y = median),data = nonmem.summary,colour = "blue",size = 1,linetype = "dashed")
		plotobj <- plotobj + scale_y_log10("Infliximab Concentration (mg/L)\n",breaks = c(0.1,1,10,100),labels = c(0.1,1,10,100))
		plotobj <- plotobj + scale_x_continuous("\nTime (days)",breaks = c(0,14,42,seq(from = 98,to = 546,by = 56)),labels = c(0,14,42,seq(from = 98,to = 546,by = 56)))
		plotobj

# ------------------------------------------------------------------------------
# Prepare NONMEM simulation output ready for Bayesian estimation of individual PK parameters
# Remove columns not necessary for input
	nonmem.input.data <- nonmem.data[,-c(5:11,20:22)]
# Rename the ID column ready for NONMEM
	names(nonmem.input.data)[1] <- "CID"
# Re-add EVID and CMT columns
	nonmem.input.data$EVID <- 0
	nonmem.input.data$EVID[nonmem.input.data$RATE == -2] <- 1
	nonmem.input.data$CMT <- 1
# Assign sample times
	sample.times <- c(98,322,546)
	nonmem.input.data$DV[!nonmem.input.data$TIME %in% sample.times] <- "."	# If not a sample time, remove the DV value
	nonmem.input.data$DV[nonmem.input.data$TIME %in% sample.times & nonmem.input.data$RATE == -2] <- "."
	nonmem.input.data$MDV <- 1	# Make a MDV column
	nonmem.input.data$MDV[nonmem.input.data$TIME %in% sample.times & !c(nonmem.input.data$RATE == -2)] <- 0
	write.csv(nonmem.input.data,file = paste0(work.dir,"nonmem_bayes_input.csv"),na = ".",quote = F,row.names = F) # Write to .csv file
	
	
		