#R Script for simulating infliximab concentrations for the population
#Clinical simulation = "clinical" scenario where dosing is guided by measured trough concentration
#If the trough = 1.5 mg/L and the target trough is 3 mg/L, then double the dose
#Source the popoulation file
#----------------------------------------------------------------------------------
#####################
##_SECOND INTERVAL_##
#####################
#Simulate some error on the previous concentration - this is what dosing decisions with be based on
set.seed(111111)
ERR1 <- 1+rnorm(n,mean = 0,sd = ERRPRO)
sim.data1$DV[sim.data1$TIME == sample.time1] <- sim.data1$CONC[sim.data1$TIME == sample.time1]*ERR1
sample.sim.data1 <- sim.data1[sim.data1$TIME == sample.time1,]
input.sim.data.clinical2 <- merge(pop.data[pop.data$TIME >= sample.time1,],sample.sim.data1[c("SIM","ID","A1","A2","AUT","DV")],by = c("SIM","ID"),all = T) #Source time-points for the second interval and merge with DV sample

#If the previous concentration was within range, continue with previous dose
#If the previous SAMPLED concentration was below 3 or greater than or equal to 5 mg/L, then increase the dose by x-fold
input.sim.data.clinical2$DOSE[input.sim.data.clinical2$DV < trough.target | input.sim.data.clinical2$DV >= trough.upper] <- input.sim.data.clinical2$DOSE[input.sim.data.clinical2$DV < trough.target | input.sim.data.clinical2$DV >= trough.upper]*(trough.target/input.sim.data.clinical2$DV[input.sim.data.clinical2$DV < trough.target | input.sim.data.clinical2$DV >= trough.upper])
input.sim.data.clinical2$DOSE[!c(input.sim.data.clinical2$TIME %in% TIMEi)] <- 0

#Simulate concentrations for the second sampling interval
sim.data.clinical2 <- ddply(input.sim.data.clinical2, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)
sim.data.clinical2 <- sim.data.clinical2[sim.data.clinical2$TIME %in% TIME2all[TIME2all < max(TIME2all)],]	#Remove the last time-point because some of my other code was lazy and this makes the plot look ugly
#This is just a time-point marking when the infusion administered just after sampling ended

#----------------------------------------------------------------------------------
####################
##_THIRD INTERVAL_##
####################
#Simulate some error on the previous concentration - this is what dosing decisions with be based on
set.seed(222222)
ERR2 <- 1+rnorm(n,mean = 0,sd = ERRPRO)
sim.data.clinical2$DV[sim.data.clinical2$TIME == sample.time2] <- sim.data.clinical2$CONC[sim.data.clinical2$TIME == sample.time2]*ERR2
sample.sim.data2 <- sim.data.clinical2[sim.data.clinical2$TIME == sample.time2,]
input.sim.data.clinical3 <- merge(pop.data[pop.data$TIME >= sample.time2,],sample.sim.data2[c("SIM","ID","A1","A2","AUT","DV")],by = c("SIM","ID"),all = T) #Source time-points for the third interval and merge with DV sample

#If the previous concentration was within range, continue with previous dose
#If the previous SAMPLED concentration was below 3 or greater than or equal to 5 mg/L, then increase the dose by x-fold
input.sim.data.clinical3$DOSE[input.sim.data.clinical3$DV < trough.target | input.sim.data.clinical3$DV >= trough.upper] <- input.sim.data.clinical3$DOSE[input.sim.data.clinical3$DV < trough.target | input.sim.data.clinical3$DV >= trough.upper]*(trough.target/input.sim.data.clinical3$DV[input.sim.data.clinical3$DV < trough.target | input.sim.data.clinical3$DV >= trough.upper])
input.sim.data.clinical3$DOSE[!c(input.sim.data.clinical3$TIME %in% TIMEi)] <- 0

#Simulate concentrations for the third sampling interval
sim.data.clinical3 <- ddply(input.sim.data.clinical3, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)
sim.data.clinical3 <- sim.data.clinical3[sim.data.clinical3$TIME %in% TIME3all[TIME3all < max(TIME3all)],]	#Remove the last time-point because some of my other code was lazy and this makes the plot look ugly
#This is just a time-point marking when the infusion administered just after sampling ended

#----------------------------------------------------------------------------------
#####################
##_FOURTH INTERVAL_##
#####################
set.seed(333333)
ERR3 <- 1+rnorm(n,mean = 0,sd = ERRPRO)
sim.data.clinical3$DV[sim.data.clinical3$TIME == sample.time3] <- sim.data.clinical3$CONC[sim.data.clinical3$TIME == sample.time3]*ERR3
sample.sim.data3 <- sim.data.clinical3[sim.data.clinical3$TIME == sample.time3,]
input.sim.data.clinical4 <- merge(pop.data[pop.data$TIME >= sample.time3,],sample.sim.data3[c("SIM","ID","A1","A2","AUT","DV")],by = c("SIM","ID"),all = T) #Source time-points for the third interval and merge with DV sample

#If the previous concentration was within range, continue with previous dose
#If the previous SAMPLED concentration was below 4 or greater than or equal to 5 mg/L, then increase the dose by x-fold
input.sim.data.clinical4$DOSE[input.sim.data.clinical4$DV < trough.target | input.sim.data.clinical4$DV >= trough.upper] <- input.sim.data.clinical4$DOSE[input.sim.data.clinical4$DV < trough.target | input.sim.data.clinical4$DV >= trough.upper]*(trough.target/input.sim.data.clinical4$DV[input.sim.data.clinical4$DV < trough.target | input.sim.data.clinical4$DV >= trough.upper])
input.sim.data.clinical4$DOSE[!c(input.sim.data.clinical4$TIME %in% TIMEi)] <- 0

#Simulate concentrations for the third sampling interval
sim.data.clinical4 <- ddply(input.sim.data.clinical4, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)
sim.data.clinical4 <- sim.data.clinical4[sim.data.clinical4$TIME %in% TIME4all,]
sim.data.clinical4$DV <- NA

#----------------------------------------------------------------------------------
#Combine the sim.data.clinicalX data frames together
sim.data.clinical <- rbind(sim.data1[sim.data1$TIME < max(TIME1),],
													sim.data.clinical2[sim.data.clinical2$TIME < max(TIME2),],
													sim.data.clinical3[sim.data.clinical3$TIME < max(TIME3),],
													sim.data.clinical4)
													
#Add DV column back into sim.data.clinical for later reference
sim.data.clinical$DV[sim.data.clinical$TIME == sample.time1] <- sim.data.clinical$CONC[sim.data.clinical$TIME == sample.time1]*ERR1
sim.data.clinical$DV[sim.data.clinical$TIME == sample.time2] <- sim.data.clinical$CONC[sim.data.clinical$TIME == sample.time2]*ERR2
sim.data.clinical$DV[sim.data.clinical$TIME == sample.time3] <- sim.data.clinical$CONC[sim.data.clinical$TIME == sample.time3]*ERR3

#Collect the cumulative time spent under the target trough concentration at each sampling time for each simulated individual
AUT.collect.sim.data.clinical <- ddply(sim.data.clinical, .(SIM,ID), AUT.collect.function)
AUT.collect.sim.data.clinical$METHOD <- "sim.clinical"

#----------------------------------------------------------------------------------
##################
##_SAVE_RESULTS_##
##################
#Save "sim.data.clinical" to a .csv file
filename3 <- "sim.data.clinical.csv"
sim.data.clinical <- sim.data.clinical[with(sim.data.clinical, order(sim.data.clinical$SIM,sim.data.clinical$ID,sim.data.clinical$TIME)), ]
write.csv(sim.data.clinical,file = filename3,na = ".",quote = F,row.names = F)

#Remove objects from workspace to free up some memory
#Only keep AUT.collect because data keeps being added to it
rm(input.sim.data.clinical2,sim.data.clinical2,input.sim.data.clinical3,sim.data.clinical3,input.sim.data.clinical4,sim.data.clinical4,sim.data.clinical)
