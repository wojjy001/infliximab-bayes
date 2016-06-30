#R Script for simulating infliximab concentrations for the population
#Optimised simulation = doses are optimised using maximum likelihood estimation
#At the time of optimising dose, it assumed that covariate and random effects are constant for the next sampling interval
#----------------------------------------------------------------------------------
#####################
##_SECOND INTERVAL_##
#####################
#Create a data frame of the dose in sampling interval 1
dose.sampling1 <- sim.data1[c("SIM","ID","TIME","DOSE")]
dose.sampling1 <- dose.sampling1[dose.sampling1$TIME == dose.time1,]
#Create a data frame of sample and final conditions for sampling interval 1
sim.sampling1 <- sim.data1[c("SIM","ID","TIME","ADA","ALB","A1","A2","CONC","AUT")]
sim.sampling1 <- sim.sampling1[sim.sampling1$TIME == sample.time1,]
#Create a data frame of individual parameters for sampling interval 1
ind.par1 <- sim.data1[c("SIM","ID","TIME","ETA1","ETA2","ETA3","ETA4")]
ind.par1 <- ind.par1[ind.par1$TIME == sample.time1,]
ind.par1 <- ind.par1[-3]

#Create a data frame for dose optimisation
interval <- "second"
status <- "simulated"
input.optim.data1 <- create.input.optim.data.function(sim.sampling1,ind.par1,TIME2all[TIME2all < max(TIME2all)])
optim.dose1 <- ddply(input.optim.data1, .(SIM,ID), optimise.dose.function.cmpf, .parallel = TRUE)
optim.dose1$INT <- 2

#Simulate the second interval using the optimised doses
input.sim.optim.data2 <- create.input.sim.optim.function(dose.sampling1,input.optim.data1,sample.time1,optim.dose1,sim.sampling1,TIME2all[TIME2all < max(TIME2all)])
sim.data.optim2 <- ddply(input.sim.optim.data2, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)

#----------------------------------------------------------------------------------
#####################
##_THIRD INTERVAL_##
#####################
#Create a data frame of the dose in sampling interval 2
dose.sampling2 <- sim.data.optim2[c("SIM","ID","TIME","DOSE")]
dose.sampling2 <- dose.sampling2[dose.sampling2$TIME == dose.time2,]
#Create a data frame of sample and final conditions for sampling interval 2
sim.sampling2 <- sim.data.optim2[c("SIM","ID","TIME","ADA","ALB","A1","A2","CONC","AUT")]
sim.sampling2 <- sim.sampling2[sim.sampling2$TIME == sample.time2,]
#Create a data frame of individual parameters for sampling interval 2
ind.par2 <- sim.data.optim2[c("SIM","ID","TIME","ETA1","ETA2","ETA3","ETA4")]
ind.par2 <- ind.par2[ind.par2$TIME == sample.time2,]
ind.par2 <- ind.par2[-3]

#Create a data frame for dose optimisation
interval <- "third"
status <- "simulated"
input.optim.data2 <- create.input.optim.data.function(sim.sampling2,ind.par2,TIME3all[TIME3all < max(TIME3all)])
optim.dose2 <- ddply(input.optim.data2, .(SIM,ID), optimise.dose.function.cmpf, .parallel = TRUE)
optim.dose2$INT <- 3

#Simulate the third interval using the optimised doses
input.sim.optim.data3 <- create.input.sim.optim.function(dose.sampling2,input.optim.data2,sample.time2,optim.dose2,sim.sampling2,TIME3all[TIME3all < max(TIME3all)])
sim.data.optim3 <- ddply(input.sim.optim.data3, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)

#----------------------------------------------------------------------------------
#####################
##_FOURTH INTERVAL_##
#####################
#Create a data frame of the dose in sampling interval 3
dose.sampling3 <- sim.data.optim3[c("SIM","ID","TIME","DOSE")]
dose.sampling3 <- dose.sampling3[dose.sampling3$TIME == dose.time3,]
#Create a data frame of sample and final conditions for sampling interval 3
sim.sampling3 <- sim.data.optim3[c("SIM","ID","TIME","ADA","ALB","A1","A2","CONC","AUT")]
sim.sampling3 <- sim.sampling3[sim.sampling3$TIME == sample.time3,]
#Create a data frame of individual parameters for sampling interval 3
ind.par3 <- sim.data.optim3[c("SIM","ID","TIME","ETA1","ETA2","ETA3","ETA4")]
ind.par3 <- ind.par3[ind.par3$TIME == sample.time3,]
ind.par3 <- ind.par3[-3]

#Create a data frame for dose optimisation
interval <- "fourth"
status <- "simulated"
input.optim.data3 <- create.input.optim.data.function(sim.sampling3,ind.par3,TIME4all)
optim.dose3 <- ddply(input.optim.data3, .(SIM,ID), optimise.dose.function.cmpf, .parallel = TRUE)
optim.dose3$INT <- 4

#Simulate the fourth interval using the optimised doses
input.sim.optim.data4 <- create.input.sim.optim.function(dose.sampling3,input.optim.data3,sample.time3,optim.dose3,sim.sampling3,TIME4all)
sim.data.optim4 <- ddply(input.sim.optim.data4, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)

#----------------------------------------------------------------------------------
#Combine the sim.data.optimX data frames together
sim.data.optim1 <- sim.data1[-ncol(sim.data1)]
sim.data.optim <- rbind(sim.data.optim1[sim.data.optim1$TIME < max(TIME1),],
													sim.data.optim2[sim.data.optim2$TIME < max(TIME2),],
													sim.data.optim3[sim.data.optim3$TIME < max(TIME3),],
													sim.data.optim4)												

#Collect the cumulative time spent under the target trough concentration at each sampling time for each simulated individual
AUT.collect.sim.data.optim <- ddply(sim.data.optim, .(SIM,ID), AUT.collect.function)
AUT.collect.sim.data.optim$METHOD <- "sim.optim"

#----------------------------------------------------------------------------------
##################
##_SAVE_RESULTS_##
##################
#Save "sim.data.optim" to a .csv file
filename4 <- "sim.data.optim.csv"
sim.data.optim <- sim.data.optim[with(sim.data.optim, order(sim.data.optim$SIM,sim.data.optim$ID,sim.data.optim$TIME)), ]
write.csv(sim.data.optim,file = filename4,na = ".",quote = F,row.names = F)

#----------------------------------------------------------------------------------
#Combine AUT.collect results into a single data frame and save as .csv file
AUT.collect <- rbind(AUT.collect.sim.data.label,
										AUT.collect.sim.data.clinical,
										AUT.collect.sim.data.optim)
AUT.collect <- AUT.collect[with(AUT.collect, order(AUT.collect$SIM,AUT.collect$ID)), ]
filename <- "AUT.collect.csv"
write.csv(AUT.collect,file = filename,na = ".",quote = F,row.names = F)
