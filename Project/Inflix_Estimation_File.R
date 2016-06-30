#R Script for Bayesian forecasting infliximab concentrations for an individual
#----------------------------------------------------------------------------------
########################################
##_SAMPLE, FIT AND SIM FIRST INTERVAL_##
########################################
#Sample the first concentration on Day 98
#Sample ADA status and ALB at time of sampling
#Assume ADA status was constant throughout the sampling interval
#Assume ALB change linearly with respect to time throughout the sampling interval
#Covariate information can be found in "cov.time.data" or "pop.data"
#Estimate a set of individual random effect parameters that best describe the individual given different time-weight methods and levels of covariate information

sim.bayes.optim.data1 <- sim.data1[sim.data1$SIM != 0,]
#Make a data frame that will serve as "initial estimates" for the Bayesian function
bayes.result0 <- ID.data2
bayes.result0$ETA1 <- 0	#For the first interval, make everyone's ETA values = 0
bayes.result0$ETA2 <- 0	#i.e., assuming that are all typical population individuals
bayes.result0$ETA3 <- 0
bayes.result0$ETA4 <- 0
#Create a data frame of the dose in sampling interval 1
dose.sampling1 <- sim.bayes.optim.data1[c("SIM","ID","TIME","DOSE")]
dose.sampling1 <- dose.sampling1[dose.sampling1$TIME == dose.time1,]
#Create a data frame of compartment and AUT initial conditions for sampling interval 1
init.sampling1 <- sim.bayes.optim.data1[c("SIM","ID","TIME","A1","A2","AUT")]
init.sampling1 <- init.sampling1[init.sampling1$TIME == dose.time1,]
#Create a data frame of sample and final conditions for sampling interval 1
final.sampling1 <- sim.bayes.optim.data1[c("SIM","ID","TIME","CONC","A1","A2","AUT")]
final.sampling1 <- final.sampling1[final.sampling1$TIME == sample.time1,]
#Add error to sampling concentration at the time of sample.time1
set.seed(111111)
final.sampling1$DV <- final.sampling1$CONC*(1+rnorm(n,mean = 0,sd = ERRPRO))

#Create a data frame fit for Bayesian estimation of the first sampling interval
bTIMEseq1 <- TIME1all[TIME1all < max(TIME1all)]
input.bayes.data1 <- create.input.bayes.data.function(covariate,dose.sampling1,init.sampling1,final.sampling1,dose.time1,sample.time1,bayes.result0,bTIMEseq1,method)
interval <- "first"
bayes.result1 <- ddply(input.bayes.data1, .(SIM,ID), bayesian.function.cmpf, .parallel = TRUE)
bayes.result1$INT <- 1

#Simulate a concentration time-profile for each individual based on fitted individual parameters and "known" covariate values
input.sim.bayes.data1 <- create.input.bayesfit.data.function(input.bayes.data1,bayes.result1)
sim.data.bayes1 <- ddply(input.sim.bayes.data1, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)	#Data frame of simulated concentrations based on the individuals Bayesian estimated parameters and "known" covariate values

#Key data frames!
#bayes.result1 contains individual fitted random effect parameters
#sim.data.bayes1 contains simulated concentrations based on individual fitted random effect parameters

#----------------------------------------------------------------------------------
#######################################
##_OPTIMISE DOSE FOR SECOND INTERVAL_##
#######################################
#Based on "known" covariate values
#Based on Bayesian estimated fitted parameters
#Start from conditions from Bayesian simulation
bayesian.sampling1 <- sim.data.bayes1[c("SIM","ID","TIME","ADA","ALB","A1","A2","CONC","AUT")]
bayesian.sampling1 <- bayesian.sampling1[bayesian.sampling1$TIME == sample.time1,]

#Create a data frame for dose optimisation
interval <- "second"
status <- "Bayesian estimated"
input.optim.data1 <- create.input.optim.data.function(bayesian.sampling1,bayes.result1,TIME2all[TIME2all < max(TIME2all)])
optim.dose1 <- ddply(input.optim.data1, .(SIM,ID), optimise.dose.function.cmpf, .parallel = TRUE)
optim.dose1$INT <- 2

#Simulate concentrations for the "TRUE" individual given the "optimised" dose
#However, if the Bayesian fitting said that the individual's trough concentration was in the target range, then continue with the previous dose
input.sim.bayes.optim.data2 <- create.input.sim.optim.function(dose.sampling1,sim.data.bayes1,sample.time1,optim.dose1,final.sampling1,TIME2all[TIME2all < max(TIME2all)])
sim.bayes.optim.data2 <- ddply(input.sim.bayes.optim.data2, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)

#----------------------------------------------------------------------------------
#########################################
##_SAMPLE, FIT AND SIM SECOND INTERVAL_##
#########################################
#Sample the second concentration on Day 210
#Sample ADA status and ALB at time of sampling
#Assume ADA status was constant throughout the sampling interval
#Assume ALB change linearly with respect to time throughout the sampling interval
#Covariate information can be found in "cov.time.data" or "pop.data"
#Estimate a set of individual random effect parameters that best describe the individual given different time-weight methods and levels of covariate information

#Create a data frame of the dose in sampling interval 2
dose.sampling2 <- sim.bayes.optim.data2[c("SIM","ID","TIME","DOSE")]
dose.sampling2 <- dose.sampling2[dose.sampling2$TIME == dose.time2,]
#Create a data frame of compartment and AUT initial conditions for sampling interval 2
init.sampling2 <- sim.bayes.optim.data2[c("SIM","ID","TIME","A1","A2","AUT")]
init.sampling2 <- init.sampling2[init.sampling2$TIME == dose.time2,]
#Create a data frame of sample and final conditions for sampling interval 2
final.sampling2 <- sim.bayes.optim.data2[c("SIM","ID","TIME","CONC","A1","A2","AUT")]
final.sampling2 <- final.sampling2[final.sampling2$TIME == sample.time2,]
#Add error to sampling concentration at the time of sample.time2
set.seed(222222)
final.sampling2$DV <- final.sampling2$CONC*(1+rnorm(n,mean = 0,sd = ERRPRO))

#Create a data frame fit for Bayesian estimation of the first sampling interval
bTIMEseq2 <- TIME2all[TIME2all < max(TIME2all)]
input.bayes.data1 <- create.input.bayes.data.function(covariate,dose.sampling1,init.sampling1,final.sampling1,dose.time1,sample.time1,bayes.result1,bTIMEseq1,method)
input.bayes.data2 <- create.input.bayes.data.function(covariate,dose.sampling2,init.sampling2,final.sampling2,dose.time2,sample.time2,bayes.result1,bTIMEseq2,method)
input.bayes.data12 <- rbind(input.bayes.data1,input.bayes.data2)	#Bind the previous input.bayes.data with current one
input.bayes.data12$DOSE[input.bayes.data12$TIME == sample.time1 & is.na(input.bayes.data12$DV) == F] <- 0	#Reset the dose at time of first sample - the subsequent dose has a line below of it's own
#Estimating individual random effect parameters for the individual
interval <- "second"
bayes.result2 <- ddply(input.bayes.data12, .(SIM,ID), bayesian.function.cmpf, .parallel = TRUE)
bayes.result2$INT <- 2

#Simulate a concentration time-profile for each individual based on fitted individual parameters and "known" covariate values
input.sim.bayes.data2 <- create.input.bayesfit.data.function(input.bayes.data12,bayes.result2)
sim.data.bayes2 <- ddply(input.sim.bayes.data2, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)	#Data frame of simulated concentrations based on the individuals Bayesian estimated parameters and "known" covariate values

#Key data frames!
#bayes.result2 contains individual fitted random effect parameters
#sim.data.bayes2 contains simulated concentrations based on individual fitted random effect parameters

#----------------------------------------------------------------------------------
######################################
##_OPTIMISE DOSE FOR THIRD INTERVAL_##
######################################
#Based on "known" covariate values
#Based on Bayesian estimated fitted parameters
#Start from conditions from Bayesian simulation
bayesian.sampling2 <- sim.data.bayes2[c("SIM","ID","TIME","ADA","ALB","A1","A2","CONC","AUT")]
bayesian.sampling2 <- bayesian.sampling2[bayesian.sampling2$TIME == sample.time2,]

#Create a data frame for dose optimisation
interval <- "third"
status <- "Bayesian estimated"
input.optim.data2 <- create.input.optim.data.function(bayesian.sampling2,bayes.result2,TIME3all[TIME3all < max(TIME3all)])
optim.dose2 <- ddply(input.optim.data2, .(SIM,ID), optimise.dose.function.cmpf, .parallel = TRUE)
optim.dose2$INT <- 3

#Simulate concentrations for the "TRUE" individual given the "optimised" dose
#However, if the Bayesian fitting said that the individual's trough concentration was in the target range, then continue with the previous dose
input.sim.bayes.optim.data3 <- create.input.sim.optim.function(dose.sampling2,sim.data.bayes2,sample.time2,optim.dose2,final.sampling2,TIME3all[TIME3all < max(TIME3all)])
sim.bayes.optim.data3 <- ddply(input.sim.bayes.optim.data3, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)

#----------------------------------------------------------------------------------
#########################################
##_SAMPLE, FIT AND SIM THIRD INTERVAL_##
#########################################
#Sample the third concentration on Day 320
#Sample ADA status and ALB at time of sampling
#Assume ADA status was constant throughout the sampling interval
#Assume ALB change linearly with respect to time throughout the sampling interval
#Covariate information can be found in "cov.time.data" or "pop.data"
#Estimate a set of individual random effect parameters that best describe the individual given different time-weight methods and levels of covariate information

#Create a data frame of the dose in sampling interval 3
dose.sampling3 <- sim.bayes.optim.data3[c("SIM","ID","TIME","DOSE")]
dose.sampling3 <- dose.sampling3[dose.sampling3$TIME == dose.time3,]
#Create a data frame of compartment and AUT initial conditions for sampling interval 3
init.sampling3 <- sim.bayes.optim.data3[c("SIM","ID","TIME","A1","A2","AUT")]
init.sampling3 <- init.sampling3[init.sampling3$TIME == dose.time3,]
#Create a data frame of sample and final conditions for sampling interval 3
final.sampling3 <- sim.bayes.optim.data3[c("SIM","ID","TIME","CONC","A1","A2","AUT")]
final.sampling3 <- final.sampling3[final.sampling3$TIME == sample.time3,]
#Add error to sampling concentration at the time of sample.time3
set.seed(333333)
final.sampling3$DV <- final.sampling3$CONC*(1+rnorm(n,mean = 0,sd = ERRPRO))

#Create a data frame fit for Bayesian estimation of the first sampling interval
bTIMEseq3 <- TIME3all[TIME3all < max(TIME3all)]
input.bayes.data1 <- create.input.bayes.data.function(covariate,dose.sampling1,init.sampling1,final.sampling1,dose.time1,sample.time1,bayes.result2,bTIMEseq1,method)
input.bayes.data2 <- create.input.bayes.data.function(covariate,dose.sampling2,init.sampling2,final.sampling2,dose.time2,sample.time2,bayes.result2,bTIMEseq2,method)
input.bayes.data3 <- create.input.bayes.data.function(covariate,dose.sampling3,init.sampling3,final.sampling3,dose.time3,sample.time3,bayes.result2,bTIMEseq3,method)
input.bayes.data123 <- rbind(input.bayes.data1,input.bayes.data2,input.bayes.data3)	#Bind the previous input.bayes.data with current one
input.bayes.data123$DOSE[input.bayes.data123$TIME == sample.time1 | input.bayes.data123$TIME == sample.time2 & is.na(input.bayes.data123$DV) == F] <- 0	#Reset the dose at time of first sample - the subsequent dose has a line below of it's own
#Estimating individual random effect parameters for the individual
interval <- "third"
bayes.result3 <- ddply(input.bayes.data123, .(SIM,ID), bayesian.function.cmpf, .parallel = TRUE)
bayes.result3$INT <- 3

#Simulate a concentration time-profile for each individual based on fitted individual parameters and "known" covariate values
input.sim.bayes.data3 <- create.input.bayesfit.data.function(input.bayes.data123,bayes.result3)
sim.data.bayes3 <- ddply(input.sim.bayes.data3, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)	#Data frame of simulated concentrations based on the individuals Bayesian estimated parameters and "known" covariate values

#Key data frames!
#bayes.result3 contains individual fitted random effect parameters
#sim.data.bayes3 contains simulated concentrations based on individual fitted random effect parameters

#----------------------------------------------------------------------------------
#######################################
##_OPTIMISE DOSE FOR FOURTH INTERVAL_##
#######################################
#Based on "known" covariate values
#Based on Bayesian estimated fitted parameters
#Start from conditions from Bayesian simulation
bayesian.sampling3 <- sim.data.bayes3[c("SIM","ID","TIME","ADA","ALB","A1","A2","CONC","AUT")]
bayesian.sampling3 <- bayesian.sampling3[bayesian.sampling3$TIME == sample.time3,]

#Create a data frame for dose optimisation
interval <- "fourth"
status <- "Bayesian estimated"
input.optim.data3 <- create.input.optim.data.function(bayesian.sampling3,bayes.result3,TIME4all)
optim.dose3 <- ddply(input.optim.data3, .(SIM,ID), optimise.dose.function.cmpf, .parallel = TRUE)
optim.dose3$INT <- 4

#Simulate concentrations for the "TRUE" individual given the "optimised" dose
#However, if the Bayesian fitting said that the individual's trough concentration was in the target range, then continue with the previous dose
input.sim.bayes.optim.data4 <- create.input.sim.optim.function(dose.sampling3,sim.data.bayes3,sample.time3,optim.dose3,final.sampling3,TIME4all)
sim.bayes.optim.data4 <- ddply(input.sim.bayes.optim.data4, .(SIM,ID), simulate.conc.cmpf, .parallel = TRUE)

#Create a data frame of the dose in sampling interval 4
dose.sampling4 <- sim.bayes.optim.data4[c("SIM","ID","TIME","DOSE")]
dose.sampling4 <- dose.sampling4[dose.sampling4$TIME == dose.time4,]
#Create a data frame of compartment and AUT initial conditions for sampling interval 4
init.sampling4 <- sim.bayes.optim.data4[c("SIM","ID","TIME","A1","A2","AUT")]
init.sampling4 <- init.sampling4[init.sampling4$TIME == dose.time4,]

#----------------------------------------------------------------------------------
#Combine the sim.bayes.optim.dataX data frames together
sim.bayes.optim.data1 <- sim.data1[-ncol(sim.data1)]
sim.bayes.optim.data1 <- sim.bayes.optim.data1[sim.bayes.optim.data1$SIM != 0,]
sim.bayes.optim.data <- rbind(sim.bayes.optim.data1[sim.bayes.optim.data1$TIME < max(TIME1),],
													sim.bayes.optim.data2[sim.bayes.optim.data2$TIME < max(TIME2),],
													sim.bayes.optim.data3[sim.bayes.optim.data3$TIME < max(TIME3),],
													sim.bayes.optim.data4)

# #When I don't run the clinical simulation file, do this instead
# sim.data1 <- sim.data1[sim.data1$SIM != 0,]													
# sim.bayes.optim.data <- rbind(sim.data1[sim.data1$TIME < max(TIME1),],
													# sim.bayes.optim.data2[sim.bayes.optim.data2$TIME < max(TIME2),],
													# sim.bayes.optim.data3[sim.bayes.optim.data3$TIME < max(TIME3),],
													# sim.bayes.optim.data4)

#Collect the cumulative time spent under the target trough concentration at each sampling time for each simulated individual
AUT.collect.bayes.data.optim <- ddply(sim.bayes.optim.data, .(SIM,ID), AUT.collect.function)
AUT.collect.bayes.data.optim$METHOD <- paste(method,covariate,sep = "")


