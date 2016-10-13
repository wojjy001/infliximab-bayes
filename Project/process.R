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
  n <- 9  # Number of seed individuals that were simulated
  nsim <- 10  # Number of simulations of seed individuals per set
  nset <- length(file.list)  # Number of simulation sets in the directory
  set.seq <- 1:nset  # Sequence of "set" numbers
  input.list <- data.frame(file.list,set.seq)

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
  levels(all.data$STUDYf) <- c("Label","Clinical","Clinical TDM","Bayes - 5mg/kg Init","Bayes - 10mg/kg Init")

# Give each individual a unique ID number (uID)
  uID <- sort(c(rep(1:9000,times = length(all.data$IPRE)/9000)))
  all.data$uID <- uID

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
  ind.summary.data <- ddply(all.data, .(STUDY,set.seq,SIM,ID,uID), ind.summary.function)

# Line (median) and ribbon (prediction intervals)
  # Bin time
    ind.summary.data$TIMEBIN <- ind.summary.data$time
    ind.summary.data$TIMEBIN[ind.summary.data$time > 98 & ind.summary.data$time <= 126] <- 98
    ind.summary.data$TIMEBIN[ind.summary.data$time > 126 & ind.summary.data$time <= 182] <- 154
    ind.summary.data$TIMEBIN[ind.summary.data$time > 182 & ind.summary.data$time <= 238] <- 210
    ind.summary.data$TIMEBIN[ind.summary.data$time > 238 & ind.summary.data$time <= 294] <- 266
    ind.summary.data$TIMEBIN[ind.summary.data$time > 294 & ind.summary.data$time <= 350] <- 322
    ind.summary.data$TIMEBIN[ind.summary.data$time > 350 & ind.summary.data$time <= 406] <- 378
    ind.summary.data$TIMEBIN[ind.summary.data$time > 406 & ind.summary.data$time <= 462] <- 434
    ind.summary.data$TIMEBIN[ind.summary.data$time > 462 & ind.summary.data$time <= 518] <- 490
    ind.summary.data$TIMEBIN[ind.summary.data$time > 518 & ind.summary.data$time <= 600] <- 546

  # Assign descriptions to ID
    ind.summary.data$IPRE <- as.numeric(ind.summary.data$IPRE)
    ind.summary.data$WT <- as.numeric(ind.summary.data$WT)
    ind.summary.data$ALB <- as.numeric(ind.summary.data$ALB)
    ind.summary.data$IDf <- as.factor(ind.summary.data$ID)
    levels(ind.summary.data$IDf) <- c("WT 40, ALB 2.5","WT 40, ALB 3","WT 40, ALB 3.5","WT 70, ALB 2.5","WT 70, ALB 3","WT 70, ALB 3.5","WT 100, ALB 2.5","WT 100, ALB 3","WT 100, ALB 3.5")
  # Assign descriptions to STUDY
    ind.summary.data$STUDYf <- as.factor(ind.summary.data$STUDY)
    levels(ind.summary.data$STUDYf) <- c("Label","Clinical","Clinical TDM","Bayes - 5mg/kg Init","Bayes - 10mg/kg Init")

## Plot concentration-time for all studies
# Calculate
  STUDY.summary.data <- ddply(ind.summary.data, .(STUDY,STUDYf,TIMEBIN), function(ind.summary.data) summary.function(ind.summary.data$IPRE))
  ID.STUDY.summary.data <- ddply(ind.summary.data, .(ID,IDf,STUDY,STUDYf,TIMEBIN), function(ind.summary.data) summary.function(ind.summary.data$IPRE))

# Plot by STUDY
  plotobj1 <- NULL
  plotobj1 <- ggplot(STUDY.summary.data)
  # plotobj1 <- plotobj1 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.95lo,ymax = stat.95hi,fill = STUDYf),alpha = 0.3)
  plotobj1 <- plotobj1 + geom_line(aes(x = time,y = IPRE,colour = STUDYf,group = uID),data = ind.summary.data,alpha = 0.05)
  plotobj1 <- plotobj1 + geom_line(aes(x = TIMEBIN,y = median),data = STUDY.summary.data)
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 3),linetype = "dashed")
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 5),linetype = "dashed")
  plotobj1 <- plotobj1 + scale_y_log10("Trough Infliximab Concentrations (mg/L)\n",breaks = c(0.001,0.01,0.1,1,10,100),labels = c(0.001,0.01,0.1,1,10,100),lim = c(0.0001,1000))
  plotobj1 <- plotobj1 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
  plotobj1 <- plotobj1 + facet_wrap(~STUDYf,ncol = 5)
  plotobj1 <- plotobj1 + theme(legend.position = "none")
  plotobj1

  ggsave(plot = plotobj1,filename = paste0(project.dir,"trough_time_noCI.png"),units = "cm",width = 30,height = 10)

# Plot - Each ID separately by STUDY
  ID.list <- 1:9
  plot.function <- function(ID.list) {
    plotobj2 <- NULL
    plotobj2 <- ggplot(ID.STUDY.summary.data[ID.STUDY.summary.data$ID == ID.list,])
    plotobj2 <- plotobj2 + ggtitle(ID.STUDY.summary.data$IDf[ID.STUDY.summary.data$ID == ID.list])
    plotobj2 <- plotobj2 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.95lo,ymax = stat.95hi,fill = STUDYf),alpha = 0.3)
    plotobj2 <- plotobj2 + geom_line(aes(x = TIMEBIN,y = median,colour = STUDYf))
    plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 3),linetype = "dashed")
    plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 5),linetype = "dashed")
    plotobj2 <- plotobj2 + scale_y_log10("Trough Infliximab Concentrations (mg/L)\n",breaks = c(0.001,0.01,0.1,1,10,100),labels = c(0.001,0.01,0.1,1,10,100))
    plotobj2 <- plotobj2 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
    plotobj2 <- plotobj2 + facet_wrap(~STUDYf,ncol = 5)
    plotobj2 <- plotobj2 + theme(legend.position = "none")
    plotobj2

    ID.label <- ID.STUDY.summary.data$IDf[ID.STUDY.summary.data$ID == ID.list][1]
    ggsave(plot = plotobj2,filename = paste0(project.dir,ID.label,"_trough_time.png"),units = "cm",width = 30,height = 10)
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
  plotobj3 <- plotobj3 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.95lo,ymax = stat.95hi,fill = STUDYf),alpha = 0.3)
  plotobj3 <- plotobj3 + geom_line(aes(x = TIMEBIN,y = median,colour = STUDYf))
  plotobj3 <- plotobj3 + scale_y_continuous("Weight at Trough Times (kg)\n")
  plotobj3 <- plotobj3 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
  plotobj3 <- plotobj3 + facet_wrap(~STUDYf,ncol = 5)
  plotobj3 <- plotobj3 + theme(legend.position = "none")
  plotobj3

  ggsave(plot = plotobj3,filename = paste0(project.dir,"WT_time.png"),units = "cm",width = 30,height = 10)

# Plot - Each ID separately by STUDY
  wt.plot.function <- function(ID.list) {
    plotobj4 <- NULL
    plotobj4 <- ggplot(ID.STUDY.wt.summary.data[ID.STUDY.wt.summary.data$ID == ID.list,])
    plotobj4 <- plotobj4 + ggtitle(ID.STUDY.wt.summary.data$IDf[ID.STUDY.wt.summary.data$ID == ID.list])
    plotobj4 <- plotobj4 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.95lo,ymax = stat.95hi,fill = STUDYf),alpha = 0.3)
    plotobj4 <- plotobj4 + geom_line(aes(x = TIMEBIN,y = median,colour = STUDYf))
    plotobj4 <- plotobj4 + scale_y_continuous("Weight at Trough Times (kg)\n")
    plotobj4 <- plotobj4 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
    plotobj4 <- plotobj4 + facet_wrap(~STUDYf,ncol = 5)
    plotobj4 <- plotobj4 + theme(legend.position = "none")
    plotobj4

    ID.label <- ID.STUDY.wt.summary.data$IDf[ID.STUDY.wt.summary.data$ID == ID.list][1]
    ggsave(plot = plotobj4,filename = paste0(project.dir,ID.label,"_WT_time.png"),units = "cm",width = 30,height = 10)
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
  plotobj5 <- plotobj5 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.95lo,ymax = stat.95hi,fill = STUDYf),alpha = 0.3)
  plotobj5 <- plotobj5 + geom_line(aes(x = TIMEBIN,y = median,colour = STUDYf))
  plotobj5 <- plotobj5 + scale_y_continuous("Albumin at Trough Times (U/L)\n")
  plotobj5 <- plotobj5 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
  plotobj5 <- plotobj5 + facet_wrap(~STUDYf,ncol = 5)
  plotobj5 <- plotobj5 + theme(legend.position = "none")
  plotobj5

  ggsave(plot = plotobj5,filename = paste0(project.dir,"ALB_time.png"),units = "cm",width = 30,height = 10)

# Plot - Each ID separately by STUDY
  alb.plot.function <- function(ID.list) {
    plotobj6 <- NULL
    plotobj6 <- ggplot(ID.STUDY.alb.summary.data[ID.STUDY.alb.summary.data$ID == ID.list,])
    plotobj6 <- plotobj6 + ggtitle(ID.STUDY.alb.summary.data$IDf[ID.STUDY.alb.summary.data$ID == ID.list])
    plotobj6 <- plotobj6 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.95lo,ymax = stat.95hi,fill = STUDYf),alpha = 0.3)
    plotobj6 <- plotobj6 + geom_line(aes(x = TIMEBIN,y = median,colour = STUDYf))
    plotobj6 <- plotobj6 + scale_y_continuous("Albumin at Trough Times (U/L)\n")
    plotobj6 <- plotobj6 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
    plotobj6 <- plotobj6 + facet_wrap(~STUDYf,ncol = 5)
    plotobj6 <- plotobj6 + theme(legend.position = "none")
    plotobj6

    ID.label <- ID.STUDY.alb.summary.data$IDf[ID.STUDY.alb.summary.data$ID == ID.list][1]
    ggsave(plot = plotobj6,filename = paste0(project.dir,ID.label,"_ALB_time.png"),units = "cm",width = 30,height = 10)
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
  plotobj7 <- plotobj7 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.95lo,ymax = stat.95hi,fill = STUDYf),alpha = 0.3)
  plotobj7 <- plotobj7 + geom_line(aes(x = TIMEBIN,y = median,colour = STUDYf))
  plotobj7 <- plotobj7 + geom_hline(aes(yintercept = 0.1),linetype = "dashed")
  plotobj7 <- plotobj7 + scale_y_continuous("Proportion of time under target trough\n")
  plotobj7 <- plotobj7 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
  plotobj7 <- plotobj7 + facet_wrap(~STUDYf,ncol = 5)
  plotobj7 <- plotobj7 + theme(legend.position = "none")
  plotobj7

  ggsave(plot = plotobj7,filename = paste0(project.dir,"ptut_time.png"),units = "cm",width = 30,height = 10)

# Plot - Each ID separately by STUDY
  ptut.plot.function <- function(ID.list) {
    plotobj8 <- NULL
    plotobj8 <- ggplot(ID.STUDY.ptut.summary.data[ID.STUDY.ptut.summary.data$ID == ID.list,])
    plotobj8 <- plotobj8 + ggtitle(ID.STUDY.ptut.summary.data$IDf[ID.STUDY.ptut.summary.data$ID == ID.list])
    plotobj8 <- plotobj8 + geom_ribbon(aes(x = TIMEBIN,ymin = stat.95lo,ymax = stat.95hi,fill = STUDYf),alpha = 0.3)
    plotobj8 <- plotobj8 + geom_line(aes(x = TIMEBIN,y = median,colour = STUDYf))
    plotobj8 <- plotobj8 + geom_hline(aes(yintercept = 0.1),linetype = "dashed")
    plotobj8 <- plotobj8 + scale_y_continuous("Proportion of time under target trough\n")
    plotobj8 <- plotobj8 + scale_x_continuous("\nTime (days)",breaks = c(0,98,210,322,434,546),labels = c(0,98,210,322,434,546))
    plotobj8 <- plotobj8 + facet_wrap(~STUDYf,ncol = 5)
    plotobj8 <- plotobj8 + theme(legend.position = "none")
    plotobj8

    ID.label <- ID.STUDY.ptut.summary.data$IDf[ID.STUDY.ptut.summary.data$ID == ID.list][1]
    ggsave(plot = plotobj8,filename = paste0(project.dir,ID.label,"_ptut_time.png"),units = "cm",width = 30,height = 10)
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

# Change in total body weight, kg
  last.data <- all.data[all.data$time == 546,]
  last.data$cWT <- last.data$WTCOV-last.data$BASE_WT
  wt.summary <- ddply(last.data, .(STUDY), function(last.data) summary.function(last.data$cWT))

# Change in albumin, U/L
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
  ada.on.proportion <- ddply(ADA.onset.summary[ADA.onset.summary$time.on != 600,], .(STUDY,STUDYf), function(ADA.onset.summary) length(ADA.onset.summary$time.on)/9000)
  ada.revert.proportion <- ddply(ADA.onset.summary[ADA.onset.summary$time.on != 600 & ADA.onset.summary$time.off != 600,], .(STUDY), function(ADA.onset.summary) length(ADA.onset.summary$time.off)/9000)

## Study Numerical Summaries by Seed
# Proportion of time under target trough at study conclusion, days (pTUT)
  ID.pTUT.summary <- ddply(all.data[all.data$time == 546,], .(STUDY,STUDYf,ID,IDf), function(all.data) summary.function(all.data$pTUT))

# Maintenance dose, mg/kg
  ID.mg.kg.summary <- ddply(main.data, .(STUDY,STUDYf,ID,IDf), function(main.data) summary.function(main.data$mg.kg))

# Maintenance dose frequency, days
  ID.int.summary <- ddply(main.data, .(STUDY,STUDYf,ID,IDf), function(main.data) summary.function(main.data$int))

# Change in total body weight, kg
  ID.wt.summary <- ddply(last.data, .(STUDY,STUDYf,ID,IDf), function(last.data) summary.function(last.data$cWT))

# Change in albumin, U/L
  ID.alb.summary <- ddply(last.data, .(STUDY,STUDYf,ID,IDf), function(last.data) summary.function(last.data$cALB))

# Proportion that develop ADA
  ID.ada.on.summary <- ddply(ADA.onset.summary[ADA.onset.summary$time.on != 600,], .(STUDY,STUDYf,ID,IDf), function(ADA.onset.summary) summary.function(ADA.onset.summary$time.on))
  ID.ada.on.proportion <- ddply(ADA.onset.summary[ADA.onset.summary$time.on != 600,], .(STUDY,STUDYf,ID,IDf), function(ADA.onset.summary) length(ADA.onset.summary$time.on)/1000)
  ID.ada.revert.proportion <- ddply(ADA.onset.summary[ADA.onset.summary$time.on != 600 & ADA.onset.summary$time.off != 600,], .(STUDY,STUDYf,ID,IDf), function(ADA.onset.summary) length(ADA.onset.summary$time.off)/1000)

## Study graphical summaries by seed
# Proportion of time below target trough
  plotobj9 <- NULL
  plotobj9 <- ggplot(ID.pTUT.summary)
  plotobj9 <- plotobj9 + geom_point(aes(x = IDf,y = median,colour = STUDYf),size = 4,alpha = 0.3)
  plotobj9 <- plotobj9 + geom_point(aes(x = IDf,y = median,colour = STUDYf),size = 4,shape = 1)
  plotobj9 <- plotobj9 + scale_x_discrete("\nBaseline seed")
  plotobj9 <- plotobj9 + scale_y_continuous("Median proportion of time below target trough\n")
  plotobj9 <- plotobj9 + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
  plotobj9

  ggsave(plot = plotobj9,filename = paste0(project.dir,"pTUT_546_Seed.png"),units = "cm",width = 20,height = 20)

# Change in weight
  plotobj10 <- NULL
  plotobj10 <- ggplot(ID.wt.summary)
  plotobj10 <- plotobj10 + geom_point(aes(x = IDf,y = median,colour = STUDYf),size = 4,alpha = 0.3)
  plotobj10 <- plotobj10 + geom_point(aes(x = IDf,y = median,colour = STUDYf),size = 4,shape = 1)
  plotobj10 <- plotobj10 + scale_x_discrete("\nBaseline seed")
  plotobj10 <- plotobj10 + scale_y_continuous("Median change in weight (kg)\n")
  plotobj10 <- plotobj10 + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
  plotobj10

  ggsave(plot = plotobj10,filename = paste0(project.dir,"cWT_Seed.png"),units = "cm",width = 20,height = 20)

# Change in albumin
  plotobj11 <- NULL
  plotobj11 <- ggplot(ID.alb.summary)
  plotobj11 <- plotobj11 + geom_point(aes(x = IDf,y = median,colour = STUDYf),size = 4,alpha = 0.3)
  plotobj11 <- plotobj11 + geom_point(aes(x = IDf,y = median,colour = STUDYf),size = 4,shape = 1)
  plotobj11 <- plotobj11 + scale_x_discrete("\nBaseline seed")
  plotobj11 <- plotobj11 + scale_y_continuous("Median change in albumin (U/L)\n")
  plotobj11 <- plotobj11 + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
  plotobj11

  ggsave(plot = plotobj11,filename = paste0(project.dir,"cALB_Seed.png"),units = "cm",width = 20,height = 20)

# Proportion with ADA onset and proportion that revert back
  ID.ada.on.proportion$ON <- 0
  ID.ada.revert.proportion$ON <- 1
  ID.ada <- rbind(ID.ada.on.proportion,ID.ada.revert.proportion)
  ID.ada$ON <- as.factor(ID.ada$ON)
  levels(ID.ada$ON) <- c("Onset of ADA","Reverted ADA status")

  plotobj12 <- NULL
  plotobj12 <- ggplot(ID.ada)
  plotobj12 <- plotobj12 + geom_point(aes(x = IDf,y = V1,colour = STUDYf),size = 4,alpha = 0.3)
  plotobj12 <- plotobj12 + geom_point(aes(x = IDf,y = V1,colour = STUDYf),size = 4,shape = 1)
  plotobj12 <- plotobj12 + scale_x_discrete("\nBaseline seed")
  plotobj12 <- plotobj12 + scale_y_continuous("Proportion of individuals\n")
  plotobj12 <- plotobj12 + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
  plotobj12 <- plotobj12 + facet_wrap(~ON)
  plotobj12

  ggsave(plot = plotobj12,filename = paste0(project.dir,"ADA_Seed.png"),units = "cm",width = 40,height = 20)
