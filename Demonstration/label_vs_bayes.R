# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))
# Load package libraries
  library(ggplot2)  # Plotting
  library(grid)  # Plotting
  library(plyr)  # ddply function

# Custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 10))
	theme_bw2 <- theme_update(
		plot.title = element_text(face = "bold",hjust = 0)
	)

# Pre-specified functions
  CI95lo <- function(x) quantile(x,probs = 0.025)  # Calculate 2.5th percentile
  CI95hi <- function(x) quantile(x,probs = 0.975)  # Calculate 97.5th percentile

# Summary function for calculating median and prediction intervals
  summary.function <- function(x) {
    median <- median(x)
    stat.95lo <- quantile(x,probs = 0.025)
    stat.95hi <- quantile(x,probs = 0.975)
    result <- c(median,stat.95lo,stat.95hi)
    names(result)[c(1,2,3)] <- c("median","stat.95lo","stat.95hi")
    result
  }

# Set trough target objects
  trough.target <- 3
  trough.upper <- 5

# ------------------------------------------------------------------------------
# Set working directory
  # project.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/Moved-Infliximab-Output/" # Mac directory
  project.dir <- "E:/Wojciechowski/Moved-Infliximab-Output/"  # Server directory

# Read in simulation output
  file.list <- list.files(path = project.dir,pattern = "SUCCESS") # List of successful simulation sets
  nset <- length(file.list)  # Number of simulation sets in the directory
  set.seq <- 1:nset  # Sequence of "set" numbers
  input.list <- data.frame(file.list,set.seq)
	rand.set <- sample(set.seq,1)
	# rand.set <- 46
	input.list <- input.list[input.list$set.seq == rand.set,]

# Read in the simulation data
  read.data.function <- function(input.list) {
    # Label data
      label.data <- read.csv(file = paste0(project.dir,input.list$file.list,"/label_simulation.csv"))
      label.data$STUDY <- 1
    # Clinical data
      clinical.data <- read.csv(file = paste0(project.dir,input.list$file.list,"/clinical_simulation.csv"))
      clinical.data$STUDY <- 2
    # Clinical TDM data
      clinical.TDM.data <- read.csv(file = paste0(project.dir,input.list$file.list,"/clinical_TDM_simulation.csv"))
      clinical.TDM.data$STUDY <- 3
    # Optimise.bayes1 data (5 mg/kg initiation)
      optimise.bayes.data1 <- read.csv(file = paste0(project.dir,input.list$file.list,"/optimise_bayes_data1.csv"))
      optimise.bayes.data1$STUDY <- 4
    # Optimise.bayes2.data (10 mg/kg initiation)
      optimise.bayes.data2 <- read.csv(file = paste0(project.dir,input.list$file.list,"/optimise_bayes_data2.csv"))
      optimise.bayes.data2$STUDY <- 5
    # Bind data from the set
      set.data <- rbind(label.data,clinical.data,clinical.TDM.data,optimise.bayes.data1,optimise.bayes.data2)
  }
  all.data <- ddply(input.list, .(set.seq), read.data.function)
  all.data <- all.data[all.data$time <= 546,] # Remove NA rows
  all.data$IPRE <- as.numeric(levels(all.data$IPRE))[all.data$IPRE]

  all.data$IDf <- as.factor(all.data$ID)
  levels(all.data$IDf) <- c("WT 40, ALB 2.5","WT 40, ALB 3","WT 40, ALB 3.5","WT 70, ALB 2.5","WT 70, ALB 3","WT 70, ALB 3.5","WT 100, ALB 2.5","WT 100, ALB 3","WT 100, ALB 3.5")
# Assign descriptions to STUDY
  all.data$STUDYf <- as.factor(all.data$STUDY)
  levels(all.data$STUDYf) <- c("Label","Clinical","Clinical TDM","Bayesian Adaptive","Bayes - 10mg/kg Init")

# There are some really small negative numbers in the dataset
  all.data$IPRE[all.data$IPRE < 0] <- 1e-8

# For each individual calculate when doses were given and how much
  ind.summary.function <- function(all.data) {
    if (!c(all.data$STUDY[1] %in% c(1,3))) {
      ind.summary.data <- data.frame(
        time = all.data$time[all.data$amt != 0],
        amt = all.data$amt[all.data$amt != 0],
        int = c(0,diff(all.data$time[all.data$amt != 0])),
        IPRE = all.data$IPRE[all.data$amt != 0],
        ADA = all.data$ADA[all.data$amt != 0],
        ALB = all.data$ALBCOV[all.data$amt != 0],
        WT = all.data$WTCOV[all.data$amt != 0],
        pTUT = all.data$pTUT[all.data$amt != 0],
        "mg/kg" = all.data$amt[all.data$amt != 0]/all.data$WTCOV[all.data$amt != 0]
      )
    } else {
      ind.summary.data <- data.frame(
        time = all.data$time[all.data$amt != 0 | all.data$time == 546],
        amt = all.data$amt[all.data$amt != 0 | all.data$time == 546],
        int = c(0,diff(all.data$time[all.data$amt != 0 | all.data$time == 546])),
        IPRE = all.data$IPRE[all.data$amt != 0 | all.data$time == 546],
        ADA = all.data$ADA[all.data$amt != 0 | all.data$time == 546],
        ALB = all.data$ALBCOV[all.data$amt != 0 | all.data$time == 546],
        WT = all.data$WTCOV[all.data$amt != 0 | all.data$time == 546],
        pTUT = all.data$pTUT[all.data$amt != 0 | all.data$time == 546],
        "mg/kg" = all.data$amt[all.data$amt != 0 | all.data$time == 546]/all.data$WTCOV[all.data$amt != 0 | all.data$time == 546]
      )
    }
  }
  ind.summary.data <- ddply(all.data, .(STUDY,set.seq,SIM,ID), ind.summary.function)

# ------------------------------------------------------------------------------
# Pull out a random individual with ID == 1 (i.e., baseline weight 40 kg and baseline albumin 2.5 U/L)
	rand.sim <- sample(1:max(all.data$SIM),1)
	# rand.sim <- 7
	# rand.id <- sample(1:9,1)
	rand.id <- 9
	rand.data <- all.data[all.data$ID == rand.id & all.data$SIM == rand.sim,]

	set.name <- input.list$file.list[input.list$set.seq == rand.set]

	filename <- paste0(project.dir,"Plots/",set.name,"_SIM",rand.sim,"_ID",rand.id,"_label_versus_Bayes.png")
	png(filename,width = 3200,height = 2400,pointsize = 10,res = 300)

	#4 ggplot2 graphs in a grid layout
	vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(4,4)))

	# Plot IPRE over time for STUDY == 1 and STUDY == 4
	 	plotobj1 <- NULL
		plotobj1 <- ggplot(rand.data[rand.data$STUDY == 1 | rand.data$STUDY == 4,])
		plotobj1 <- plotobj1 + ggtitle("(a)\n")
		plotobj1 <- plotobj1 + geom_line(aes(x = time,y = IPRE,colour = STUDYf))
		plotobj1 <- plotobj1 + geom_hline(aes(yintercept = trough.target),linetype = "dashed")
		plotobj1 <- plotobj1 + geom_hline(aes(yintercept = trough.upper),linetype = "dashed")
		plotobj1 <- plotobj1 + scale_y_log10("Infliximab Concentration (mg/L)\n")
		plotobj1 <- plotobj1 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
		plotobj1 <- plotobj1 + facet_wrap(~STUDYf)
		plotobj1 <- plotobj1 + theme(legend.position = "none")
		print(plotobj1,vp = vplayout(1:2,1:2))

	# Plot weight over time
		plotobj2 <- NULL
		plotobj2 <- ggplot(rand.data[rand.data$STUDY == 1 | rand.data$STUDY == 4,])
		plotobj2 <- plotobj2 + ggtitle("(b)\n")
		plotobj2 <- plotobj2 + geom_line(aes(x = time,y = WTCOV,colour = STUDYf))
		plotobj2 <- plotobj2 + scale_y_continuous("Total Body Weight (kg)\n")
		plotobj2 <- plotobj2 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
		plotobj2 <- plotobj2 + facet_wrap(~STUDYf)
		plotobj2 <- plotobj2 + theme(legend.position = "none")
		print(plotobj2,vp = vplayout(1:2,3:4))

	# Plot albumin over time
		plotobj3 <- NULL
		plotobj3 <- ggplot(rand.data[rand.data$STUDY == 1 | rand.data$STUDY == 4,])
		plotobj3 <- plotobj3 + ggtitle("(c)\n")
		plotobj3 <- plotobj3 + geom_line(aes(x = time,y = ALBCOV,colour = STUDYf))
		plotobj3 <- plotobj3 + scale_y_continuous("Albumin (U/L)\n")
		plotobj3 <- plotobj3 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
		plotobj3 <- plotobj3 + facet_wrap(~STUDYf)
		plotobj3 <- plotobj3 + theme(legend.position = "none")
		print(plotobj3,vp = vplayout(3:4,1:2))

	# Plot proportion of time under target trough over time
		plotobj4 <- NULL
		plotobj4 <- ggplot(rand.data[rand.data$STUDY == 1 | rand.data$STUDY == 4,])
		plotobj4 <- plotobj4 + ggtitle("(d)\n")
		plotobj4 <- plotobj4 + geom_line(aes(x = time,y = pTUT,colour = STUDYf))
		plotobj4 <- plotobj4 + scale_y_continuous("Proportion of time under target trough\n")
		plotobj4 <- plotobj4 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
		plotobj4 <- plotobj4 + facet_wrap(~STUDYf)
		plotobj4 <- plotobj4 + theme(legend.position = "none")
		print(plotobj4,vp = vplayout(3:4,3:4))

	dev.off()

	rand.summary.data <- ind.summary.data[ind.summary.data$ID == rand.id & ind.summary.data$SIM == rand.sim & ind.summary.data$STUDY %in% c(1,4),]
