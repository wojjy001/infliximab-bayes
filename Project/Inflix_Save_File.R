#R Script for saving the results
print(paste("Scenario (",method," ",covariate,") has started.",sep = "")) 	#Print a "initiation" message
source(paste(work.dir,"Pop4_Inflix_Estimation_File.R",sep = ""))	#Source and run the estimation file
#----------------------------------------------------------------------------------
##################
##_SAVE_RESULTS_##
##################
#Save method specific results to a method specific folder
method.output.dir <- paste(sim.output.dir,method,covariate,sep = "")
dir.create(file.path(method.output.dir),showWarnings = FALSE)
setwd(file.path(method.output.dir))

#----------------------------------------------------------------------------------
#Combine dose.sampling data frames together and write to .csv file
dose.sampling <- rbind(dose.sampling1,dose.sampling2,dose.sampling3,dose.sampling4)
dose.sampling <- dose.sampling[with(dose.sampling, order(dose.sampling$SIM,dose.sampling$ID,dose.sampling$TIME)), ]
filename5 <- "dose.sampling.csv"
write.csv(dose.sampling,file = filename5,na = ".",quote = F,row.names = F)

#Combine init.sampling data frames together and write to .csv file
init.sampling <- rbind(init.sampling1,init.sampling2,init.sampling3,init.sampling4)
init.sampling <- init.sampling[with(init.sampling, order(init.sampling$SIM,init.sampling$ID,init.sampling$TIME)), ]
filename6 <- "init.sampling.csv"
write.csv(init.sampling,file = filename6,na = ".",quote = F,row.names = F)

#Combine final.sampling data frames together and write. .csv file
final.sampling <- rbind(final.sampling1,final.sampling2,final.sampling3)
final.sampling <- final.sampling[with(final.sampling, order(final.sampling$SIM,final.sampling$ID,final.sampling$TIME)), ]
filename7 <- "final.sampling.csv"
write.csv(final.sampling,file = filename7,na = ".",quote = F,row.names = F)

#Combine bayesian.sampling data frames together and write. .csv file
bayesian.sampling <- rbind(bayesian.sampling1,bayesian.sampling2,bayesian.sampling3)
bayesian.sampling <- bayesian.sampling[with(bayesian.sampling, order(bayesian.sampling$SIM,bayesian.sampling$ID,bayesian.sampling$TIME)), ]
filename8 <- "bayesian.sampling.csv"
write.csv(bayesian.sampling,file = filename8,na = ".",quote = F,row.names = F)

#----------------------------------------------------------------------------------
#Write the sequential Bayes fit data frames to .csv files
#Fitted concentrations for the first sampling interval
fit.data.bayes1 <- sim.data.bayes1
fit.data.bayes1 <- fit.data.bayes1[with(fit.data.bayes1, order(fit.data.bayes1$SIM,fit.data.bayes1$ID,fit.data.bayes1$TIME)), ]
filename9 <- "fit.data.bayes1.csv"
write.csv(fit.data.bayes1,file = filename9,na = ".",quote = F,row.names = F)

#Fitted concentrations for the second sampling interval
fit.data.bayes2 <- sim.data.bayes2
fit.data.bayes2 <- fit.data.bayes2[with(fit.data.bayes2, order(fit.data.bayes2$SIM,fit.data.bayes2$ID,fit.data.bayes2$TIME)), ]
filename10 <- "fit.data.bayes2.csv"
write.csv(fit.data.bayes2,file = filename10,na = ".",quote = F,row.names = F)

#Fitted concentrations for the third sampling interval
fit.data.bayes3 <- sim.data.bayes3
fit.data.bayes3 <- fit.data.bayes3[with(fit.data.bayes3, order(fit.data.bayes3$SIM,fit.data.bayes3$ID,fit.data.bayes3$TIME)), ]
filename11 <- "fit.data.bayes3.csv"
write.csv(fit.data.bayes3,file = filename11,na = ".",quote = F,row.names = F)

#----------------------------------------------------------------------------------
#Write the Bayesian results to a .csv file (individual parameter estimates)
bayes.result <- rbind(bayes.result1,bayes.result2,bayes.result3)
bayes.result <- bayes.result[with(bayes.result, order(bayes.result$SIM,bayes.result$ID,bayes.result$INT)), ]
filename12 <- "bayes.result.csv"
write.csv(bayes.result,file = filename12,na = ".",quote = F,row.names = F)

#Write the dose optimisation results to a .csv file
optim.dose <- rbind(optim.dose1,optim.dose2,optim.dose3)
optim.dose <- optim.dose[with(optim.dose, order(optim.dose$SIM,optim.dose$ID,optim.dose$INT)), ]
filename13 <- "optim.dose.csv"
write.csv(optim.dose,file = filename13,na = ".",quote = F,row.names = F)

#----------------------------------------------------------------------------------
#Write the "TRUE" simulated concentrations as a result of Bayesian fitting individual parameters and dose optimisation to a .csv file
filename14 <- "sim.bayes.optim.data.csv"
sim.bayes.optim.data <- sim.bayes.optim.data[with(sim.bayes.optim.data, order(sim.bayes.optim.data$SIM,sim.bayes.optim.data$ID,sim.bayes.optim.data$TIME)), ]
write.csv(sim.bayes.optim.data,file = filename14,na = ".",quote = F,row.names = F)

#----------------------------------------------------------------------------------
#Combine AUT.collect results into a single data frame and save as .csv file
setwd(sim.output.dir)
AUT.collect <- rbind(AUT.collect,
										AUT.collect.bayes.data.optim)
AUT.collect <- AUT.collect[with(AUT.collect, order(AUT.collect$SIM,AUT.collect$ID,AUT.collect$AUT11)), ]
write.csv(AUT.collect,file = filename,na = ".",quote = F,row.names = F)

#----------------------------------------------------------------------------------
###############
##_COMPLETED_##
###############
#Remove objects from workspace to free up some memory
#Only keep AUT.collect because data keeps being added to it
#Sample, fit and sim first interval
rm(dose.sampling1,init.sampling1,final.sampling1,input.bayes.data1,bayes.result1,input.sim.bayes.data1,sim.data.bayes1)
#Optimise dose for second interval
rm(bayesian.sampling1,input.optim.data1,optim.dose1,input.sim.bayes.optim.data2)
#Sample, fit and sim second interval
rm(dose.sampling2,init.sampling2,final.sampling2,input.bayes.data2,input.bayes.data12,bayes.result2,input.sim.bayes.data2,sim.data.bayes2)
#Optimise dose for third interval
rm(bayesian.sampling2,input.optim.data2,optim.dose2,input.sim.bayes.optim.data3,sim.bayes.optim.data3)
#Sample, fit and sim third interval
rm(dose.sampling3,init.sampling3,final.sampling3,input.bayes.data3,input.bayes.data123,bayes.result3,input.sim.bayes.data3,sim.data.bayes3)
#Optimise dose for fourth interval
rm(bayesian.sampling3,input.optim.data3,optim.dose3,input.sim.bayes.optim.data4,sim.bayes.optim.data4,dose.sampling4,init.sampling4)
#Save results
rm(dose.sampling,init.sampling,final.sampling,bayesian.sampling,fit.data.bayes1,fit.data.bayes2,fit.data.bayes3,bayes.result.optim.dose,sim.bayes.optim.data)

#Print a "completed" message
print(paste("Scenario (",method," ",covariate,") has finished.",sep = ""))
setwd(sim.output.dir)
