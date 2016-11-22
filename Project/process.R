# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))
# Load package libraries
  library(ggplot2)  # Plotting
  library(grid)  # Plotting
  library(plyr)  # ddply function

# Custom ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))

# Pre-specified functions
  CI95lo <- function(x) quantile(x,probs = 0.025)  # Calculate 2.5th percentile
  CI95hi <- function(x) quantile(x,probs = 0.975)  # Calculate 97.5th percentile

# Summary function for calculating median and prediction intervals
  summary.function <- function(x) {
    # stat.95lo <- quantile(x,probs = 0.025,names = F)
		stat.90lo <- quantile(x,probs = 0.05,names = F)
		stat.80lo <- quantile(x,probs = 0.1,names = F)
		stat.60lo <- quantile(x,probs = 0.2,names = F)
		stat.50lo <- quantile(x,probs = 0.25,names = F)
		stat.40lo <- quantile(x,probs = 0.3,names = F)
		stat.20lo <- quantile(x,probs = 0.4,names = F)
    median <- median(x)
		stat.20hi <- quantile(x,probs = 0.6,names = F)
		stat.40hi <- quantile(x,probs = 0.7,names = F)
		stat.50hi <- quantile(x,probs = 0.75,names = F)
		stat.60hi <- quantile(x,probs = 0.8,names = F)
		stat.80hi <- quantile(x,probs = 0.9,names = F)
		stat.90hi <- quantile(x,probs = 0.95,names = F)
    # stat.95hi <- quantile(x,probs = 0.975,names = F)
		n <- length(x)
    result <- c(median,stat.20lo,stat.20hi,stat.40lo,stat.40hi,stat.50lo,stat.50hi,
			stat.60lo,stat.60hi,stat.80lo,stat.80hi,stat.90lo,stat.90hi,n)
		names(result) <- c("median","stat.20lo","stat.20hi","stat.40lo","stat.40hi","stat.50lo","stat.50hi",
			"stat.60lo","stat.60hi","stat.80lo","stat.80hi","stat.90lo","stat.90hi","n")
    result
  }

# Set trough target objects
  trough.target <- 3
  trough.upper <- 5

	n <- 9
	nsim <- 1000

# ------------------------------------------------------------------------------
# Set working directory
  project.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/Moved-Infliximab-Output/" # Mac directory
  # project.dir <- "E:/Wojciechowski/Moved-Infliximab-Output/"  # Server directory
	plot.dir <- paste0(project.dir,"Plots3/")

# Read in simulation output
  file.list <- list.files(path = project.dir,pattern = "SUCCESS") # List of successful simulation sets
  nset <- length(file.list)  # Number of simulation sets in the directory
  set.seq <- 1:nset  # Sequence of "set" numbers
  input.list <- data.frame(file.list,set.seq)
	input.list <- head(input.list,nsim/10)

# Read in the simulation data
  read.data.function <- function(input.list) {
    # Label data
      label.data0 <- read.csv(file = paste0(project.dir,input.list$file.list,"/time_dep_0_label_simulation.csv"))
      label.data0$STUDY <- 1
			# label.data1 <- read.csv(file = paste0(project.dir,input.list$file.list,"/time_dep_1_label_simulation.csv"))
      # label.data1$STUDY <- 5
    # Clinical data
      clinical.data0 <- read.csv(file = paste0(project.dir,input.list$file.list,"/time_dep_0_clinical_simulation.csv"))
      clinical.data0$STUDY <- 2
			# clinical.data1 <- read.csv(file = paste0(project.dir,input.list$file.list,"/time_dep_1_clinical_simulation.csv"))
      # clinical.data1$STUDY <- 6
    # Clinical TDM data
      clinical.TDM.data0 <- read.csv(file = paste0(project.dir,input.list$file.list,"/time_dep_0_clinical_TDM_simulation.csv"))
      clinical.TDM.data0$STUDY <- 3
			# clinical.TDM.data1 <- read.csv(file = paste0(project.dir,input.list$file.list,"/time_dep_1_clinical_TDM_simulation.csv"))
      # clinical.TDM.data1$STUDY <- 7
    # Optimise.bayes1 data (5 mg/kg initiation)
      optimise.bayes.data0 <- read.csv(file = paste0(project.dir,input.list$file.list,"/time_dep_0_optimise_bayes_data1.csv"))
      optimise.bayes.data0$STUDY <- 4
			# optimise.bayes.data1 <- read.csv(file = paste0(project.dir,input.list$file.list,"/time_dep_1_optimise_bayes_data1.csv"))
      # optimise.bayes.data1$STUDY <- 8
    # Bind data from the set
      set.data <- rbind(label.data0,clinical.data0,clinical.TDM.data0,optimise.bayes.data0)
			# set.data <- rbind(label.data1,clinical.data1,clinical.TDM.data1,optimise.bayes.data1)
  }
  all.data <- ddply(input.list, .(set.seq), read.data.function)
  all.data <- all.data[all.data$time <= 546,] # Remove NA rows
  # all.data$IPRE <- as.numeric(levels(all.data$IPRE))[all.data$IPRE]

  all.data$IDf <- as.factor(all.data$ID)
  levels(all.data$IDf) <- c("WT 40, ALB 2.5","WT 40, ALB 3","WT 40, ALB 3.5","WT 70, ALB 2.5","WT 70, ALB 3","WT 70, ALB 3.5","WT 100, ALB 2.5","WT 100, ALB 3","WT 100, ALB 3.5")
# Assign descriptions to STUDY
  all.data$STUDYf <- as.factor(all.data$STUDY)
  levels(all.data$STUDYf) <- c(
		# "Non-TD Label","Non-TD Clinical","Non-TD Clinical TDM","Non-TD Bayes",
	"Label","Clinical Protocol","Clinical TDM","Bayesian-Guided")

# Give each individual a unique ID number (uID)
  uID <- sort(c(rep(seq(from = 1,to = n*nsim*4,by = 1),times = length(unique(all.data$time)))))
  all.data$uID <- uID

# There are some really small negative numbers in the dataset
  all.data$IPRE[all.data$IPRE < 0] <- 1e-8

# For each individual calculate when doses were given and how much
  ind.summary.function <- function(all.data) {
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
  }
  ind.summary.data <- ddply(all.data, .(STUDY,set.seq,SIM,ID,uID), ind.summary.function)

# Line (median) and ribbon (prediction intervals)
  # Bin time
    ind.summary.data$TIMEBIN <- ind.summary.data$time
		# ind.summary.data$TIMEBIN[ind.summary.data$TIMEBIN > 98 & ind.summary.data$STUDY %in% c(2,6)] <- ceiling((ind.summary.data$TIMEBIN[ind.summary.data$TIMEBIN > 98 & ind.summary.data$STUDY %in% c(2,6)]-98)/56)*56+98
		ind.summary.data$TIMEBIN[ind.summary.data$TIMEBIN > 98 & ind.summary.data$STUDY %in% c(2,3,4,6,7,8)] <- ceiling((ind.summary.data$TIMEBIN[ind.summary.data$TIMEBIN > 98 & ind.summary.data$STUDY %in% c(2,3,4,6,7,8)]-98)/56)*56+98

		standard.times <- c(0,14,42,98,154,210,266,322,378,434,490,546)
		# subset.ind.summary.data <- ind.summary.data[ind.summary.data$time %in% standard.times,]

  # Assign descriptions to ID
    ind.summary.data$IPRE <- as.numeric(ind.summary.data$IPRE)
    ind.summary.data$WT <- as.numeric(ind.summary.data$WT)
    ind.summary.data$ALB <- as.numeric(ind.summary.data$ALB)
    ind.summary.data$IDf <- as.factor(ind.summary.data$ID)
    levels(ind.summary.data$IDf) <- c("WT 40, ALB 2.5","WT 40, ALB 3","WT 40, ALB 3.5","WT 70, ALB 2.5","WT 70, ALB 3","WT 70, ALB 3.5","WT 100, ALB 2.5","WT 100, ALB 3","WT 100, ALB 3.5")
  # Assign descriptions to STUDY
    ind.summary.data$STUDYf <- as.factor(ind.summary.data$STUDY)
    levels(ind.summary.data$STUDYf) <- c(
			# "Non-TD Label","Non-TD Clinical","Non-TD Clinical TDM","Non-TD Bayes",
		"Label","Clinical Protocol","Clinical TDM","Bayesian-Guided")

## Plot concentration-time for all studies
# Calculate
  STUDY.summary.data <- ddply(ind.summary.data, .(STUDY,STUDYf,TIMEBIN), function(ind.summary.data) summary.function(ind.summary.data$IPRE))
  ID.STUDY.summary.data <- ddply(ind.summary.data, .(ID,IDf,STUDY,STUDYf,TIMEBIN), function(ind.summary.data) summary.function(ind.summary.data$IPRE))

# Plot by STUDY
  plotobj1 <- NULL
  plotobj1 <- ggplot(STUDY.summary.data)
  plotobj1 <- plotobj1 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.90lo,ymax = stat.90hi,fill = STUDYf),alpha = 0.2)
	plotobj1 <- plotobj1 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.80lo,ymax = stat.80hi,fill = STUDYf),alpha = 0.2)
	plotobj1 <- plotobj1 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.60lo,ymax = stat.60hi,fill = STUDYf),alpha = 0.2)
	plotobj1 <- plotobj1 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.40lo,ymax = stat.40hi,fill = STUDYf),alpha = 0.2)
	plotobj1 <- plotobj1 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.20lo,ymax = stat.20hi,fill = STUDYf),alpha = 0.2)
  plotobj1 <- plotobj1 + geom_line(aes(x = TIMEBIN,y = median),data = STUDY.summary.data)
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 3),linetype = "dashed")
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 5),linetype = "dashed")
  plotobj1 <- plotobj1 + scale_y_log10("Trough Infliximab Concentrations (mg/L)\n",breaks = c(0.001,0.01,0.1,1,10,100),labels = c(0.001,0.01,0.1,1,10,100))
  plotobj1 <- plotobj1 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
  plotobj1 <- plotobj1 + facet_wrap(~STUDYf,ncol = 4)
  plotobj1 <- plotobj1 + theme(legend.position = "none")
  plotobj1

  ggsave(plot = plotobj1,filename = paste0(plot.dir,"trough_time.png"),units = "cm",width = 30,height = 10)

# Plot - Each ID separately by STUDY
  ID.list <- 1:n
  plot.function <- function(ID.list) {
    plotobj2 <- NULL
    plotobj2 <- ggplot(ID.STUDY.summary.data[ID.STUDY.summary.data$ID == ID.list,])
    plotobj2 <- plotobj2 + ggtitle(ID.STUDY.summary.data$IDf[ID.STUDY.summary.data$ID == ID.list])
		plotobj2 <- plotobj2 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.90lo,ymax = stat.90hi,fill = STUDYf),alpha = 0.2)
		plotobj2 <- plotobj2 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.80lo,ymax = stat.80hi,fill = STUDYf),alpha = 0.2)
		plotobj2 <- plotobj2 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.60lo,ymax = stat.60hi,fill = STUDYf),alpha = 0.2)
		plotobj2 <- plotobj2 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.40lo,ymax = stat.40hi,fill = STUDYf),alpha = 0.2)
		plotobj2 <- plotobj2 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.20lo,ymax = stat.20hi,fill = STUDYf),alpha = 0.2)
    plotobj2 <- plotobj2 + geom_line(aes(x = TIMEBIN,y = median))
    plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 3),linetype = "dashed")
    plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 5),linetype = "dashed")
    plotobj2 <- plotobj2 + scale_y_log10("Trough Infliximab Concentrations (mg/L)\n",breaks = c(0.001,0.01,0.1,1,10,100),labels = c(0.001,0.01,0.1,1,10,100))
    plotobj2 <- plotobj2 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
    plotobj2 <- plotobj2 + facet_wrap(~STUDYf,ncol = 4)
    plotobj2 <- plotobj2 + theme(legend.position = "none")
    plotobj2

    ID.label <- ID.STUDY.summary.data$IDf[ID.STUDY.summary.data$ID == ID.list][1]
    ggsave(plot = plotobj2,filename = paste0(plot.dir,ID.label,"_trough_time.png"),units = "cm",width = 30,height = 10)
  }
  plot.list <- lapply(ID.list,plot.function)
  plot.list

## Plot weight-time for all studies
# Calculate
  STUDY.wt.summary.data <- ddply(ind.summary.data, .(STUDY,STUDYf,TIMEBIN), function(ind.summary.data) summary.function(ind.summary.data$WT))
  ID.STUDY.wt.summary.data <- ddply(ind.summary.data, .(ID,IDf,STUDY,STUDYf,TIMEBIN), function(ind.summary.data) summary.function(ind.summary.data$WT))

# Plot by STUDY
  plotobj3 <- NULL
  plotobj3 <- ggplot(STUDY.wt.summary.data)
	plotobj3 <- plotobj3 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.90lo,ymax = stat.90hi,fill = STUDYf),alpha = 0.2)
	plotobj3 <- plotobj3 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.80lo,ymax = stat.80hi,fill = STUDYf),alpha = 0.2)
	plotobj3 <- plotobj3 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.60lo,ymax = stat.60hi,fill = STUDYf),alpha = 0.2)
	plotobj3 <- plotobj3 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.40lo,ymax = stat.40hi,fill = STUDYf),alpha = 0.2)
	plotobj3 <- plotobj3 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.20lo,ymax = stat.20hi,fill = STUDYf),alpha = 0.2)
  plotobj3 <- plotobj3 + geom_line(aes(x = TIMEBIN,y = median))
  plotobj3 <- plotobj3 + scale_y_continuous("Weight at Trough Times (kg)\n")
  plotobj3 <- plotobj3 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
  plotobj3 <- plotobj3 + facet_wrap(~STUDYf,ncol = 4)
  plotobj3 <- plotobj3 + theme(legend.position = "none")
  plotobj3

  ggsave(plot = plotobj3,filename = paste0(plot.dir,"WT_time.png"),units = "cm",width = 30,height = 10)

# Plot - Each ID separately by STUDY
  wt.plot.function <- function(ID.list) {
    plotobj4 <- NULL
    plotobj4 <- ggplot(ID.STUDY.wt.summary.data[ID.STUDY.wt.summary.data$ID == ID.list,])
    plotobj4 <- plotobj4 + ggtitle(ID.STUDY.wt.summary.data$IDf[ID.STUDY.wt.summary.data$ID == ID.list])
		plotobj4 <- plotobj4 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.90lo,ymax = stat.90hi,fill = STUDYf),alpha = 0.2)
		plotobj4 <- plotobj4 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.80lo,ymax = stat.80hi,fill = STUDYf),alpha = 0.2)
		plotobj4 <- plotobj4 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.60lo,ymax = stat.60hi,fill = STUDYf),alpha = 0.2)
		plotobj4 <- plotobj4 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.40lo,ymax = stat.40hi,fill = STUDYf),alpha = 0.2)
		plotobj4 <- plotobj4 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.20lo,ymax = stat.20hi,fill = STUDYf),alpha = 0.2)
    plotobj4 <- plotobj4 + geom_line(aes(x = TIMEBIN,y = median))
    plotobj4 <- plotobj4 + scale_y_continuous("Weight at Trough Times (kg)\n")
    plotobj4 <- plotobj4 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
    plotobj4 <- plotobj4 + facet_wrap(~STUDYf,ncol = 4)
    plotobj4 <- plotobj4 + theme(legend.position = "none")
    plotobj4

    ID.label <- ID.STUDY.wt.summary.data$IDf[ID.STUDY.wt.summary.data$ID == ID.list][1]
    ggsave(plot = plotobj4,filename = paste0(plot.dir,ID.label,"_WT_time.png"),units = "cm",width = 30,height = 10)
  }
  wt.plot.list <- lapply(ID.list,wt.plot.function)
  wt.plot.list

## Plot albumin-time for all studies
# Calculate
  STUDY.alb.summary.data <- ddply(ind.summary.data, .(STUDY,STUDYf,TIMEBIN), function(ind.summary.data) summary.function(ind.summary.data$ALB))
  ID.STUDY.alb.summary.data <- ddply(ind.summary.data, .(ID,IDf,STUDY,STUDYf,TIMEBIN), function(ind.summary.data) summary.function(ind.summary.data$ALB))

# Plot by STUDY
  plotobj5 <- NULL
  plotobj5 <- ggplot(STUDY.alb.summary.data)
	plotobj5 <- plotobj5 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.90lo,ymax = stat.90hi,fill = STUDYf),alpha = 0.2)
	plotobj5 <- plotobj5 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.80lo,ymax = stat.80hi,fill = STUDYf),alpha = 0.2)
	plotobj5 <- plotobj5 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.60lo,ymax = stat.60hi,fill = STUDYf),alpha = 0.2)
	plotobj5 <- plotobj5 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.40lo,ymax = stat.40hi,fill = STUDYf),alpha = 0.2)
	plotobj5 <- plotobj5 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.20lo,ymax = stat.20hi,fill = STUDYf),alpha = 0.2)
  plotobj5 <- plotobj5 + geom_line(aes(x = TIMEBIN,y = median))
  plotobj5 <- plotobj5 + scale_y_continuous("Albumin at Trough Times (g/dL)\n")
  plotobj5 <- plotobj5 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
  plotobj5 <- plotobj5 + facet_wrap(~STUDYf,ncol = 4)
  plotobj5 <- plotobj5 + theme(legend.position = "none")
  plotobj5

  ggsave(plot = plotobj5,filename = paste0(plot.dir,"ALB_time.png"),units = "cm",width = 30,height = 10)

# Plot - Each ID separately by STUDY
  alb.plot.function <- function(ID.list) {
    plotobj6 <- NULL
    plotobj6 <- ggplot(ID.STUDY.alb.summary.data[ID.STUDY.alb.summary.data$ID == ID.list,])
    plotobj6 <- plotobj6 + ggtitle(ID.STUDY.alb.summary.data$IDf[ID.STUDY.alb.summary.data$ID == ID.list])
		plotobj6 <- plotobj6 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.90lo,ymax = stat.90hi,fill = STUDYf),alpha = 0.2)
		plotobj6 <- plotobj6 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.80lo,ymax = stat.80hi,fill = STUDYf),alpha = 0.2)
		plotobj6 <- plotobj6 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.60lo,ymax = stat.60hi,fill = STUDYf),alpha = 0.2)
		plotobj6 <- plotobj6 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.40lo,ymax = stat.40hi,fill = STUDYf),alpha = 0.2)
		plotobj6 <- plotobj6 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.20lo,ymax = stat.20hi,fill = STUDYf),alpha = 0.2)
    plotobj6 <- plotobj6 + geom_line(aes(x = TIMEBIN,y = median))
    plotobj6 <- plotobj6 + scale_y_continuous("Albumin at Trough Times (g/dL)\n")
    plotobj6 <- plotobj6 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
    plotobj6 <- plotobj6 + facet_wrap(~STUDYf,ncol = 4)
    plotobj6 <- plotobj6 + theme(legend.position = "none")
    plotobj6

    ID.label <- ID.STUDY.alb.summary.data$IDf[ID.STUDY.alb.summary.data$ID == ID.list][1]
    ggsave(plot = plotobj6,filename = paste0(plot.dir,ID.label,"_ALB_time.png"),units = "cm",width = 30,height = 10)
  }
  alb.plot.list <- lapply(ID.list,alb.plot.function)
  alb.plot.list

## Proportion of time below target trough concentration
# Calculate
  STUDY.ptut.summary.data <- ddply(ind.summary.data, .(STUDY,STUDYf,TIMEBIN), function(ind.summary.data) summary.function(ind.summary.data$pTUT))
  ID.STUDY.ptut.summary.data <- ddply(ind.summary.data, .(ID,IDf,STUDY,STUDYf,TIMEBIN), function(ind.summary.data) summary.function(ind.summary.data$pTUT))

# Plot by STUDY
  plotobj7 <- NULL
  plotobj7 <- ggplot(STUDY.ptut.summary.data)
	plotobj7 <- plotobj7 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.90lo,ymax = stat.90hi,fill = STUDYf),alpha = 0.2)
	plotobj7 <- plotobj7 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.80lo,ymax = stat.80hi,fill = STUDYf),alpha = 0.2)
	plotobj7 <- plotobj7 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.60lo,ymax = stat.60hi,fill = STUDYf),alpha = 0.2)
	plotobj7 <- plotobj7 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.40lo,ymax = stat.40hi,fill = STUDYf),alpha = 0.2)
	plotobj7 <- plotobj7 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.20lo,ymax = stat.20hi,fill = STUDYf),alpha = 0.2)
  plotobj7 <- plotobj7 + geom_line(aes(x = TIMEBIN,y = median))
  # plotobj7 <- plotobj7 + geom_hline(aes(yintercept = 0.1),linetype = "dashed")
  plotobj7 <- plotobj7 + scale_y_continuous("Proportion of time under target trough\n")
  plotobj7 <- plotobj7 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
  plotobj7 <- plotobj7 + facet_wrap(~STUDYf,ncol = 4)
  plotobj7 <- plotobj7 + theme(legend.position = "none")
  plotobj7

  ggsave(plot = plotobj7,filename = paste0(plot.dir,"ptut_time.png"),units = "cm",width = 30,height = 10)

# Plot - Each ID separately by STUDY
  ptut.plot.function <- function(ID.list) {
    plotobj8 <- NULL
    plotobj8 <- ggplot(ID.STUDY.ptut.summary.data[ID.STUDY.ptut.summary.data$ID == ID.list,])
    plotobj8 <- plotobj8 + ggtitle(ID.STUDY.ptut.summary.data$IDf[ID.STUDY.ptut.summary.data$ID == ID.list])
		plotobj8 <- plotobj8 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.90lo,ymax = stat.90hi,fill = STUDYf),alpha = 0.2)
		plotobj8 <- plotobj8 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.80lo,ymax = stat.80hi,fill = STUDYf),alpha = 0.2)
		plotobj8 <- plotobj8 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.60lo,ymax = stat.60hi,fill = STUDYf),alpha = 0.2)
		plotobj8 <- plotobj8 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.40lo,ymax = stat.40hi,fill = STUDYf),alpha = 0.2)
		plotobj8 <- plotobj8 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.20lo,ymax = stat.20hi,fill = STUDYf),alpha = 0.2)
    plotobj8 <- plotobj8 + geom_line(aes(x = TIMEBIN,y = median))
    # plotobj8 <- plotobj8 + geom_hline(aes(yintercept = 0.1),linetype = "dashed")
    plotobj8 <- plotobj8 + scale_y_continuous("Proportion of time under target trough\n")
    plotobj8 <- plotobj8 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
    plotobj8 <- plotobj8 + facet_wrap(~STUDYf,ncol = 4)
    plotobj8 <- plotobj8 + theme(legend.position = "none")
    plotobj8

    ID.label <- ID.STUDY.ptut.summary.data$IDf[ID.STUDY.ptut.summary.data$ID == ID.list][1]
    ggsave(plot = plotobj8,filename = paste0(plot.dir,ID.label,"_ptut_time.png"),units = "cm",width = 30,height = 10)
  }
  ptut.plot.list <- lapply(ID.list,ptut.plot.function)
  ptut.plot.list

## Study Numerical Summaries
# Proportion of time under target trough at study conclusion, days (pTUT)
  pTUT.summary <- ddply(all.data[all.data$time == 546,], .(STUDY), function(all.data) summary.function(all.data$pTUT))

# Maintenance dose, mg/kg
  main.data <- ind.summary.data[ind.summary.data$time >= 98 & ind.summary.data$amt != 0,] # Subset only maintenance dosing data
  mg.kg.summary <- ddply(main.data, .(STUDY), function(main.data) summary.function(main.data$mg.kg))

# Maintenance dose frequency, days
  int.summary <- ddply(main.data, .(STUDY), function(main.data) summary.function(main.data$int))

# Maintenance dose amounts administered
	amount.summary <- ddply(main.data, .(STUDY), function(main.data) sum(main.data$amt))

# Change in total body weight, kg
  last.data <- all.data[all.data$time == 546,]
  last.data$cWT <- last.data$WTCOV-last.data$BASE_WT
  wt.summary <- ddply(last.data, .(STUDY), function(last.data) summary.function(last.data$cWT))

# Change in albumin, g/dL
  last.data$cALB <- last.data$ALBCOV-last.data$BASE_ALB
  alb.summary <- ddply(last.data, .(STUDY), function(last.data) summary.function(last.data$cALB))

# Proportion that develop ADA
  ADA.onset.function <- function(ind.summary.data) {
    ADA.onoff <- c(0,diff(ind.summary.data$ADA))
    times <- ind.summary.data$time
    if (length(ADA.onoff[ADA.onoff != 0]) != 0) time.on <- times[ADA.onoff == 1]
    if (length(ADA.onoff[ADA.onoff != 0]) == 0) time.on <- 600
    if (length(ADA.onoff[ADA.onoff == -1]) == 1) time.off <- times[ADA.onoff == -1]
    if (length(ADA.onoff[ADA.onoff == -1]) != 1) time.off <- 600
    ADA.onset.data <- data.frame(time.on,time.off)
  }
  ADA.onset.summary <- ddply(ind.summary.data, .(STUDY,STUDYf,set.seq,SIM,ID,IDf), ADA.onset.function)
  ada.on.summary <- ddply(ADA.onset.summary[ADA.onset.summary$time.on != 600,], .(STUDY,STUDYf), function(ADA.onset.summary) summary.function(ADA.onset.summary$time.on))
  ada.on.proportion <- ddply(ADA.onset.summary[ADA.onset.summary$time.on != 600,], .(STUDY,STUDYf), function(ADA.onset.summary) length(ADA.onset.summary$time.on)/(n*nsim))
  # ada.revert.proportion <- ddply(ADA.onset.summary, .(STUDY), function(ADA.onset.summary) length(ADA.onset.summary$time.off[ADA.onset.summary$time.on != 600 & ADA.onset.summary$time.off != 600])/length(ADA.onset.summary$time.on[ADA.onset.summary$time.on != 600]))
	# ada.revert.proportion

	ada.init.summary <- ddply(ADA.onset.summary[ADA.onset.summary$time.on <= 98,], .(STUDY,STUDYf), function(ADA.onset.summary) length(ADA.onset.summary$time.on))
	ada.init.proportion <- ddply(ADA.onset.summary[ADA.onset.summary$time.on <= 98,], .(STUDY,STUDYf), function(ADA.onset.summary) length(ADA.onset.summary$time.on)/(n*nsim))
	ada.main.proportion <- (ada.on.summary$n-ada.init.summary$V1)/(n*nsim)

## Study Numerical Summaries by Seed
# Proportion of time under target trough at study conclusion, days (pTUT)
  ID.pTUT.summary <- ddply(all.data[all.data$time == 546,], .(STUDY,STUDYf,ID,IDf), function(all.data) summary.function(all.data$pTUT))

# Maintenance dose, mg/kg
  ID.mg.kg.summary <- ddply(main.data, .(STUDY,STUDYf,ID,IDf), function(main.data) summary.function(main.data$mg.kg))

# Maintenance dose frequency, days
  ID.int.summary <- ddply(main.data, .(STUDY,STUDYf,ID,IDf), function(main.data) summary.function(main.data$int))

# Change in total body weight, kg
  ID.wt.summary <- ddply(last.data, .(STUDY,STUDYf,ID,IDf), function(last.data) summary.function(last.data$cWT))

# Change in albumin, g/dL
  ID.alb.summary <- ddply(last.data, .(STUDY,STUDYf,ID,IDf), function(last.data) summary.function(last.data$cALB))

# Proportion that develop ADA
  ID.ada.on.summary <- ddply(ADA.onset.summary[ADA.onset.summary$time.on != 600,], .(STUDY,STUDYf,ID,IDf), function(ADA.onset.summary) summary.function(ADA.onset.summary$time.on))
  ID.ada.on.proportion <- ddply(ADA.onset.summary[ADA.onset.summary$time.on != 600,], .(STUDY,STUDYf,ID,IDf), function(ADA.onset.summary) length(ADA.onset.summary$time.on)/nsim)

	ID.ada.init.proportion <- ddply(ADA.onset.summary[ADA.onset.summary$time.on <= 98,], .(STUDY,STUDYf,ID,IDf), function(ADA.onset.summary) length(ADA.onset.summary$time.on)/(nsim))
	ID.ada.main.proportion <- ID.ada.init.proportion
	ID.ada.main.proportion$V1 <- (ID.ada.on.summary$n/nsim)-ID.ada.init.proportion$V1

  # ID.ada.revert.proportion <- ddply(ADA.onset.summary, .(STUDY,STUDYf,ID,IDf), function(ADA.onset.summary) length(ADA.onset.summary$time.off[ADA.onset.summary$time.on != 600 & ADA.onset.summary$time.off != 600])/length(ADA.onset.summary$time.on[ADA.onset.summary$time.on != 600]))

## Study graphical summaries by seed
# Proportion of time below target trough
  plotobj9 <- NULL
  plotobj9 <- ggplot(ID.pTUT.summary)
	# plotobj9 <- plotobj9 + geom_linerange(aes(x = IDf,ymin = stat.90lo,ymax = stat.90hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj9 <- plotobj9 + geom_linerange(aes(x = IDf,ymin = stat.80lo,ymax = stat.80hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj9 <- plotobj9 + geom_linerange(aes(x = IDf,ymin = stat.60lo,ymax = stat.60hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj9 <- plotobj9 + geom_linerange(aes(x = IDf,ymin = stat.40lo,ymax = stat.40hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj9 <- plotobj9 + geom_linerange(aes(x = IDf,ymin = stat.20lo,ymax = stat.20hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	plotobj9 <- plotobj9 + geom_errorbar(aes(x = IDf,ymin = stat.50lo,ymax = stat.50hi,colour = STUDYf),position = position_dodge(width = 0.9),size = 1)
  plotobj9 <- plotobj9 + geom_point(aes(x = IDf,y = median,colour = STUDYf),position = position_dodge(width = 0.9),size = 4)
  plotobj9 <- plotobj9 + geom_point(aes(x = IDf,y = median,group = STUDYf),position = position_dodge(width = 0.9),size = 4,shape = 1)
  plotobj9 <- plotobj9 + scale_x_discrete("\nBaseline Seed")
  plotobj9 <- plotobj9 + scale_y_continuous("Proportion of time below target trough\n",breaks = seq(from = 0,to = 1,by = 0.1))
  plotobj9 <- plotobj9 + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
	plotobj9 <- plotobj9 + theme(legend.position = "none")
	# plotobj9 <- plotobj9 + facet_wrap(~WT)
  plotobj9

  ggsave(plot = plotobj9,filename = paste0(plot.dir,"pTUT_546_Seed.png"),units = "cm",width = 20,height = 20)

# Proportion of time below target trough during maintenance phase
	main.tut.function <- function(all.data) {
		init.TUT <- all.data$TUT[all.data$time == 98]
		all.data$mainTUT <- all.data$TUT
		all.data$mainTUT[all.data$time >= 98] <- all.data$TUT[all.data$time >= 98]-init.TUT
		all.data$mainpTUT <- 0
		all.data$mainpTUT[all.data$time > 98] <- all.data$mainTUT[all.data$time > 98]/(all.data$time[all.data$time > 98]-98)
		all.data
	}
	main.ind.tut.data <- ddply(all.data[all.data$time %in% standard.times,], .(uID), main.tut.function)
	main.ind.tut.summary <- ddply(main.ind.tut.data[main.ind.tut.data$time == 546,], .(STUDY,STUDYf), function(main.ind.tut.data) summary.function(main.ind.tut.data$mainpTUT))

  ID.main.pTUT.summary <- ddply(main.ind.tut.data[main.ind.tut.data$time == 546,], .(STUDY,STUDYf,ID,IDf), function(main.ind.tut.data) summary.function(main.ind.tut.data$mainpTUT))
	ID.main.pTUT.summary$ALB[ID.main.pTUT.summary$ID %in% c(1,4,7)] <- 2.5
	ID.main.pTUT.summary$ALB[ID.main.pTUT.summary$ID %in% c(2,5,8)] <- 3
	ID.main.pTUT.summary$ALB[ID.main.pTUT.summary$ID %in% c(3,6,9)] <- 3.5
	ID.main.pTUT.summary$WT[ID.main.pTUT.summary$ID %in% c(1,2,3)] <- 40
	ID.main.pTUT.summary$WT[ID.main.pTUT.summary$ID %in% c(4,5,6)] <- 70
	ID.main.pTUT.summary$WT[ID.main.pTUT.summary$ID %in% c(7,8,9)] <- 100

	ID.main.pTUT.summary$ALBf <- as.factor(ID.main.pTUT.summary$ALB)
	levels(ID.main.pTUT.summary$ALBf) <- c("ALB 2.5 g/dL","ALB 3 g/dL","ALB 3.5 g/dL")
	ID.main.pTUT.summary$WTf <- as.factor(ID.main.pTUT.summary$WT)
	levels(ID.main.pTUT.summary$WTf) <- c("WT 40 kg","WT 70 kg","WT 100 kg")

	ID.main.pTUT.summary$ADApro <- round(ID.ada.main.proportion$V1,digits = 2)

	plotobj13 <- NULL
  plotobj13 <- ggplot(ID.main.pTUT.summary)
	# plotobj13 <- plotobj13 + geom_linerange(aes(x = IDf,ymin = stat.90lo,ymax = stat.90hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj13 <- plotobj13 + geom_linerange(aes(x = IDf,ymin = stat.80lo,ymax = stat.80hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj13 <- plotobj13 + geom_linerange(aes(x = IDf,ymin = stat.60lo,ymax = stat.60hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj13 <- plotobj13 + geom_linerange(aes(x = IDf,ymin = stat.40lo,ymax = stat.40hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj13 <- plotobj13 + geom_linerange(aes(x = IDf,ymin = stat.20lo,ymax = stat.20hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	plotobj13 <- plotobj13 + geom_errorbar(aes(x = STUDYf,ymin = stat.50lo,ymax = stat.50hi,colour = STUDYf),position = position_dodge(width = 0.2),size = 1)
  plotobj13 <- plotobj13 + geom_point(aes(x = STUDYf,y = median,colour = STUDYf),position = position_dodge(width = 0.2),size = 4)
  plotobj13 <- plotobj13 + geom_point(aes(x = STUDYf,y = median,group = STUDYf),position = position_dodge(width = 0.2),size = 4,shape = 1)
	plotobj13 <- plotobj13 + geom_text(aes(x = STUDYf,y = stat.50hi+0.07,label = ADApro))
  plotobj13 <- plotobj13 + scale_x_discrete("\nStudy")
  plotobj13 <- plotobj13 + scale_y_continuous("Proportion of time below target trough during maintenance\n",breaks = seq(from = 0,to = 1,by = 0.1))
  plotobj13 <- plotobj13 + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
	# plotobj13 <- plotobj13 + theme(legend.position = "none")
	plotobj13 <- plotobj13 + facet_grid(ALBf~WTf)
  plotobj13

  ggsave(plot = plotobj13,filename = paste0(plot.dir,"mainpTUT_546_Seed.png"),units = "cm",width = 20,height = 20)

# Change in weight
  plotobj10 <- NULL
  plotobj10 <- ggplot(ID.wt.summary)
  plotobj10 <- plotobj10 + geom_point(aes(x = IDf,y = median,colour = STUDYf),size = 4,alpha = 0.5)
  plotobj10 <- plotobj10 + geom_point(aes(x = IDf,y = median,colour = STUDYf),size = 4,shape = 1)
  plotobj10 <- plotobj10 + scale_x_discrete("\nBaseline seed")
  plotobj10 <- plotobj10 + scale_y_continuous("Median change in weight (kg)\n")
  plotobj10 <- plotobj10 + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
  plotobj10

  ggsave(plot = plotobj10,filename = paste0(plot.dir,"cWT_Seed.png"),units = "cm",width = 20,height = 20)

# Change in albumin
  plotobj11 <- NULL
  plotobj11 <- ggplot(ID.alb.summary)
  plotobj11 <- plotobj11 + geom_point(aes(x = IDf,y = median,colour = STUDYf),size = 4,alpha = 0.5)
  plotobj11 <- plotobj11 + geom_point(aes(x = IDf,y = median,colour = STUDYf),size = 4,shape = 1)
  plotobj11 <- plotobj11 + scale_x_discrete("\nBaseline seed")
  plotobj11 <- plotobj11 + scale_y_continuous("Median change in albumin (g/dL)\n")
  plotobj11 <- plotobj11 + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
  plotobj11

  ggsave(plot = plotobj11,filename = paste0(plot.dir,"cALB_Seed.png"),units = "cm",width = 20,height = 20)

# Calculate time to target in the maintenance for each individual and then summarise
	time.to.target.function <- function(ind.summary.data) {
		main.data <- ind.summary.data[ind.summary.data$time >= 98,]
		for (i in 1:nrow(main.data)) {
			if (main.data$IPRE[i] >= 3) target <- T
			if (main.data$IPRE[i] < 3) target <- F
			time.to.target <- main.data$time[i]
			if (time.to.target >= 546) time.to.target <- NA
			if (target == T) break
		}
		ind.summary.data$target <- time.to.target
		ind.summary.data
	}

	time.to.target.data <- ddply(ind.summary.data, .(uID), time.to.target.function)
	time.to.target.summary <- ddply(time.to.target.data[time.to.target.data$time == 98,], .(STUDY,STUDYf), function(time.to.target.data) summary.function(na.omit(time.to.target.data$target)))
	time.to.target.ind.summary <- ddply(time.to.target.data[time.to.target.data$time == 98,], .(STUDY,STUDYf,ID,IDf), function(time.to.target.data) summary.function(na.omit(time.to.target.data$target)))
	time.to.target.ind.summary$pro <- round(time.to.target.ind.summary$n/1000,digits = 2)

	time.to.target.ind.summary$ALB[time.to.target.ind.summary$ID %in% c(1,4,7)] <- 2.5
	time.to.target.ind.summary$ALB[time.to.target.ind.summary$ID %in% c(2,5,8)] <- 3
	time.to.target.ind.summary$ALB[time.to.target.ind.summary$ID %in% c(3,6,9)] <- 3.5
	time.to.target.ind.summary$WT[time.to.target.ind.summary$ID %in% c(1,2,3)] <- 40
	time.to.target.ind.summary$WT[time.to.target.ind.summary$ID %in% c(4,5,6)] <- 70
	time.to.target.ind.summary$WT[time.to.target.ind.summary$ID %in% c(7,8,9)] <- 100

	time.to.target.ind.summary$ALBf <- as.factor(time.to.target.ind.summary$ALB)
	levels(time.to.target.ind.summary$ALBf) <- c("ALB 2.5 g/dL","ALB 3 g/dL","ALB 3.5 g/dL")
	time.to.target.ind.summary$WTf <- as.factor(time.to.target.ind.summary$WT)
	levels(time.to.target.ind.summary$WTf) <- c("WT 40 kg","WT 70 kg","WT 100 kg")

	plotobj16 <- NULL
  plotobj16 <- ggplot(time.to.target.ind.summary)
	# plotobj16 <- plotobj16 + geom_linerange(aes(x = IDf,ymin = stat.90lo,ymax = stat.90hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj16 <- plotobj16 + geom_linerange(aes(x = IDf,ymin = stat.80lo,ymax = stat.80hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj16 <- plotobj16 + geom_linerange(aes(x = IDf,ymin = stat.60lo,ymax = stat.60hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj16 <- plotobj16 + geom_linerange(aes(x = IDf,ymin = stat.40lo,ymax = stat.40hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	# plotobj16 <- plotobj16 + geom_linerange(aes(x = IDf,ymin = stat.20lo,ymax = stat.20hi,colour = STUDYf),position = position_dodge(width = 0.9),alpha = 0.2,size = 4)
	plotobj16 <- plotobj16 + geom_errorbar(aes(x = STUDYf,ymin = stat.50lo,ymax = stat.50hi,colour = STUDYf),position = position_dodge(width = 0.2),size = 1)
  plotobj16 <- plotobj16 + geom_point(aes(x = STUDYf,y = median,colour = STUDYf),position = position_dodge(width = 0.2),size = 4)
  plotobj16 <- plotobj16 + geom_point(aes(x = STUDYf,y = median,group = STUDYf),position = position_dodge(width = 0.2),size = 4,shape = 1)
	plotobj16 <- plotobj16 + geom_text(aes(x = STUDYf,y = stat.50hi+15,label = pro))
  plotobj16 <- plotobj16 + scale_x_discrete("\nStudy")
  plotobj16 <- plotobj16 + scale_y_continuous("Time to first trough target achievement (days)\n",breaks = seq(from = 98,to = 546,by = 56))
  plotobj16 <- plotobj16 + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
	# plotobj16 <- plotobj16 + theme(legend.position = "none")
	plotobj16 <- plotobj16 + facet_grid(ALBf~WTf)
  plotobj16

  ggsave(plot = plotobj16,filename = paste0(plot.dir,"mainT2T_546_Seed.png"),units = "cm",width = 20,height = 20)

# Plot 2 individual patients and their concentration-time profiles from all studies
	# 1 = 40 kg, 2.5 g/dL
	# 2 = 100 kg, 3.5 g/dL
	random.set <- sample(all.data$set.seq,1)
	low.random <- sample(all.data$SIM[all.data$ID == 1],1)
	low.random.data <- all.data[all.data$ID == 1 & all.data$SIM == low.random & all.data$set.seq == random.set,]
	low.random.data$IDf <- "Example 2: Baseline WT 40 kg, ALB 2.5 g/dL"

	high.random <- sample(all.data$SIM[all.data$ID == 9],1)
	high.random.data <- all.data[all.data$ID == 9 & all.data$SIM == high.random & all.data$set.seq == random.set,]
	high.random.data$IDf <- "Example 1: Baseline WT 100 kg, ALB 3.5 g/dL"

	random.data <- rbind(low.random.data,high.random.data)

	plotobj14 <- NULL
	plotobj14 <- ggplot(random.data)
	plotobj14 <- plotobj14 + geom_line(aes(x = time,y = IPRE,colour = STUDYf))
	plotobj14 <- plotobj14 + geom_point(aes(x = time,y = DV,colour = STUDYf),data = random.data[random.data$amt != 0,],size = 2)
	plotobj14 <- plotobj14 + geom_point(aes(x = time,y = DV),data = random.data[random.data$amt != 0,],size = 2,shape = 1)
	plotobj14 <- plotobj14 + geom_hline(aes(yintercept = 3),linetype = "dashed")
	plotobj14 <- plotobj14 + geom_hline(aes(yintercept = 5),linetype = "dashed")
	plotobj14 <- plotobj14 + scale_y_log10("Infliximab Concentration (mg/L)\n",lim = c(0.001,1000),breaks = c(0.001,0.01,0.1,1,10,100),labels = c(0.001,0.01,0.1,1,10,100))
	plotobj14 <- plotobj14 + scale_x_continuous("\nTime (days)",lim = c(98,154),breaks = seq(from = 98,to = 154,by = 7),labels = seq(from = 98,to = 154,by = 7))
	plotobj14 <- plotobj14 + theme(legend.position = "none")
	plotobj14 <- plotobj14 + facet_wrap(~IDf)
	plotobj14

  ggsave(plot = plotobj14,filename = paste0(plot.dir,"individual_concs_first_int.png"),units = "cm",width = 20,height = 10)

	plotobj15 <- NULL
	plotobj15 <- ggplot(random.data)
	plotobj15 <- plotobj15 + geom_line(aes(x = time,y = IPRE,colour = STUDYf))
	plotobj15 <- plotobj15 + geom_point(aes(x = time,y = DV,colour = STUDYf),data = random.data[random.data$amt != 0,],size = 2)
	plotobj15 <- plotobj15 + geom_point(aes(x = time,y = DV),data = random.data[random.data$amt != 0,],size = 2,shape = 1)
	plotobj15 <- plotobj15 + geom_hline(aes(yintercept = 3),linetype = "dashed")
	plotobj15 <- plotobj15 + geom_hline(aes(yintercept = 5),linetype = "dashed")
	plotobj15 <- plotobj15 + scale_y_log10("Infliximab Concentration (mg/L)\n",lim = c(0.001,1000),breaks = c(0.001,0.01,0.1,1,10,100),labels = c(0.001,0.01,0.1,1,10,100))
	plotobj15 <- plotobj15 + scale_x_continuous("\nTime (days)",lim = c(490,546),breaks = seq(from = 490,to = 546,by = 7),labels = seq(from = 490,to = 546,by = 7))
	plotobj15 <- plotobj15 + theme(legend.position = "none")
	plotobj15 <- plotobj15 + facet_wrap(~IDf)
	plotobj15

  ggsave(plot = plotobj15,filename = paste0(plot.dir,"individual_concs_last_int.png"),units = "cm",width = 20,height = 10)

# Save the workspace so I don't have to re-read all of the data frames
	save.image(file = paste0(plot.dir,"time_dependent.RData"))
