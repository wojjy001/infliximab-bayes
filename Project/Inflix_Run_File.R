#R Script for Bayesian forecasting infliximab concentrations for an individual
#----------------------------------------------------------------------------------
#First three doses are administered at 5 mg/kg at 0, 2 and 6 weeks, and individual infliximab concentrations are simulated
#The first sample (DV) is taken at Day 98 (8 weeks after the third dose at 6 weeks)
#Individual parameters are Bayes estimated using the first sample and known covariates at that time
#Based on the BAYES ESTIMATED individual parameters and covariates at the sampled time, the dose is optimised for the next two doses 8 weeks apart.
#Just like simulation, the covariates at the sampled time-point are assumed to be constant when optimising doses for the second sampling interval.
#The next two dosing intervals are simulated based on the patient's actual covariate values AND actual parameter values using the first optimised dose.
#The process is repeated for a second sample, and the next three doses 8 weeks apart.
#At 378 days (8 weeks after the 8th dose), individual parameters are bayes estimated.

#method
#NTimeWeight = No time-weighting
#Peck1.005 = Peck method, Q = 1.005
#Peck1.01 = Peck method, Q = 1.01

#covariate
#AllCov = All covariates
#NoADA = No ADA information, assume population typical, i.e., 0
#NoALB = No albumin information, assume population typical, i.e., 4
#NoCov = No covariates, assume population typical

#----------------------------------------------------------------------------------
#Remove all current objects in the workspace
rm(list = ls(all = TRUE))
work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/Pop4/"	#Working directory
#Set patient and target information for all scenarios
n <- 12	#Number of seed individuals (where each seed individual has a different set of covariate values)
nsim <- 10	#Number of simulations of the seed individuals to perform
sim.name <- paste("SIM",nsim,"_IND",n,sep = "")

#Load package libraries
library(doParallel)	#Parallel processing

#Setting up cores to run parallel processes, thus increasing speed
#Set up a cluster of cores to run the application over
cl <- makePSOCKcluster(12)
#detectCores() searches for the number of cores that the local machine has
#Contents with makePSOCKcluster brackets can be changed to a whole number if you
#want to assign an exact number

#List packages that are required to be sent to each core for the parallel process
#The foreach package always needs to be included
#This example uses the .parallel argument in ddply which calls a function that uses
#lsoda from the deSolve package
clusterEvalQ(cl, list(library(foreach),source("D:/Wojciechowski/Franklin/Pop4/Pop4_Inflix_Functions_File.R")))

#Registers the parallel backend with the foreach package (automatically loaded when doParallel is loaded)
registerDoParallel(cl)

#Create working directory for simulation (specific for number of individuals and number of simulations performed)
sim.output.dir <- paste(work.dir,sim.name,"/",sep = "")
dir.create(file.path(sim.output.dir),showWarnings = FALSE)
setwd(file.path(sim.output.dir))

#Run the simulation files and save their output
suppressPackageStartupMessages(	#Suppress package loading messages
	suppressWarnings(	#Suppress any warnings
		source(paste(work.dir,"Pop4_Inflix_Label_Simulation_File.R",sep = ""))
	)
)
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Clinical_Simulation_File.R",sep = "")))
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Optimised_Simulation_File.R",sep = "")))

#Read in AUT.collect.csv file
file.name.in <- "AUT.collect.csv"
AUT.collect <- read.csv(file.name.in,stringsAsFactors = F,na.strings = c(".","."))

#----------------------------------------------------------------------------------
#Set method and covariate information for each scenario and run sequentially
#Scenarios with No time-weighting
method <- "NTimeWeight"
covariate <- "AllCov"	#AllCov = All covariates
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))
covariate <- "NoADA"	#NoADA = No ADA information, assume population typical, i.e., 0
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))
covariate <- "NoALB"	#NoALB = No albumin information, assume population typical, i.e., 4
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))
covariate <- "NoCov"	#NoCov = No covariates, assume population typical
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))

#Scenarios using Peck, Q = 1.005
method <- "Peck1.005"
covariate <- "AllCov"	#AllCov = All covariates
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))
covariate <- "NoADA"	#NoADA = No ADA information, assume population typical, i.e., 0
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))
covariate <- "NoALB"	#NoALB = No albumin information, assume population typical, i.e., 4
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))
covariate <- "NoCov"	#NoCov = No covariates, assume population typical
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))

#Scenarios using Peck, Q = 1.01
method <- "Peck1.01"
covariate <- "AllCov"	#AllCov = All covariates
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))
covariate <- "NoADA"	#NoADA = No ADA information, assume population typical, i.e., 0
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))
covariate <- "NoALB"	#NoALB = No albumin information, assume population typical, i.e., 4
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))
covariate <- "NoCov"	#NoCov = No covariates, assume population typical
suppressWarnings(source(paste(work.dir,"Pop4_Inflix_Save_File.R",sep = "")))

#----------------------------------------------------------------------------------
#Quit R once all completed
q("no")
