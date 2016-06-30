#R Script for reading in saved output and processing the results
#----------------------------------------------------------------------------------
#Read in saved data
NTimeWeightAllCov.data <- read.csv(paste0(work.dir,"/SIM10_IND12/NTimeWeightAllCov/sim.bayes.optim.data.csv"))
plot.data <-  NTimeWeightAllCov.data[NTimeWeightAllCov.data$ID == 1 & NTimeWeightAllCov.data$SIM == 1,]

plotobj <- NULL
plotobj <- ggplot(data = NTimeWeightAllCov.data[NTimeWeightAllCov.data$ID == 1 & NTimeWeightAllCov.data$SIM == 1,])
plotobj <- plotobj + geom_line(aes(x = TIME,y = CONC),colour = "red")
plotobj <- plotobj + geom_hline(aes(yintercept = 3),linetype = "dashed")
plotobj <- plotobj + geom_hline(aes(yintercept = 5),linetype = "dashed")
plotobj <- plotobj + scale_x_continuous("\nTime (days)")
plotobj <- plotobj + scale_y_log10("Infliximab Concentration (mg/L)\n")
plotobj

POPV1 <- 3.33
#If someone's exp(ETA2) changes by 10% from baseline then what impact does that have on their V1?
BASE_ETA2 <- 0  #Assuming typical population value for V1
EXP_BASE_ETA2 <- exp(0)  #Exponential to calculate a proportional change as ETA2 is a random variable normally distributed around 0
EXP_BASE_ETA2
EXP_FINAL_ETA2 <- EXP_BASE_ETA2*1.1  #Final exp(ETA2) is 10% higher than baseline ETA2
EXP_FINAL_ETA2
FINAL_ETA2 <- log(EXP_FINAL_ETA2)

BASE_V1 <- POPV1*exp(BASE_ETA2)
FINAL_V1 <- POPV1*exp(FINAL_ETA2)
BASE_V1
FINAL_V1









