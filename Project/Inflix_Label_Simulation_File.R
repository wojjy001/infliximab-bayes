#R Script for simulating infliximab concentrations for the population
#Label simulation = every dose is 5 mg/kg regardless of trough concentrations
#Source the population file
source(paste(work.dir,"Pop4_Inflix_Population_File.R",sep = ""))
#----------------------------------------------------------------------------------
#########################
##_REMAINING INTERVALS_##
#########################
#Simulate concentrations for the entire time period where every dose is 5 mg/kg
input.sim.data.label <- pop.data
input.sim.data.label$A1 <- 0	#Central compartment initial conditions
input.sim.data.label$A2 <- 0	#Peripheral compartment initial conditions
input.sim.data.label$AUT <- 0	#Time under trough initial conditions
sim.data.label <- ddply(input.sim.data.label, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)

#Collect the cumulative time spent under the target trough concentration at each sampling time for each simulated individual
AUT.collect.sim.data.label <- ddply(sim.data.label, .(SIM,ID), AUT.collect.function)
AUT.collect.sim.data.label$METHOD <- "sim.label"

#----------------------------------------------------------------------------------
##################
##_SAVE_RESULTS_##
##################
#Save "pop.data" from "Pop4_Inflix_Population_File.R" to a .csv file
filename1 <- "pop.data.csv"
pop.data <- pop.data[with(pop.data, order(pop.data$SIM,pop.data$ID,pop.data$TIME)), ]
write.csv(pop.data,file = filename1,na = ".",quote = F,row.names = F)

#Save "sim.data.label" to a .csv file
filename2 <- "sim.data.label.csv"
sim.data.label <- sim.data.label[with(sim.data.label, order(sim.data.label$SIM,sim.data.label$ID,sim.data.label$TIME)), ]
write.csv(sim.data.label,file = filename2,na = ".",quote = F,row.names = F)

#Remove objects from workspace to free up some memory
#Only keep AUT.collect because data keeps being added to it
rm(sim.data.label)