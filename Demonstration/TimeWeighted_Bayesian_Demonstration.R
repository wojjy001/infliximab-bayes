# Create some plots that demonstrate the concept of time-weighted Bayesian estimation
# ------------------------------------------------------------------------------
# Load package libraries
	library(ggplot2) # Plotting
	library(grid) # Plotting

# Custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 12))
	theme_bw2 <- theme_update(plot.title = element_text(hjust = 0))

# Set working directory
	work.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Demonstration/"
	plotoutput.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes-plots/"

# ------------------------------------------------------------------------------
# Plot of observation's contribution to OFV over time
	ERRPRO <- 0.419	# Model's residual unexplained variability
	TIME.seq <- seq(from = 0,to = 378,by = 1) # Time sequence

	filename <- paste0(plotoutput.dir,"time_weight_demonstration_plots.png")
	png(filename,width = 900,height = 1200)
	vplayout <- function(x,y) viewport(layout.pos.row = x,layout.pos.col = y)
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(3,2)))

	NTimeWeightSD <- ERRPRO*1^(378-TIME.seq)
	Peck1005SD <- ERRPRO*1.005^(378-TIME.seq)
	Peck101SD <- ERRPRO*1.01^(378-TIME.seq)

	time.weight.sd.data <- data.frame(TIME = TIME.seq,
		SD = c(NTimeWeightSD,Peck1005SD,Peck101SD),
		Method = c(rep("NTimeWeight",times = length(TIME.seq)),rep("Peck1.005",times = length(TIME.seq)),rep("Peck1.01",times = length(TIME.seq)))
	)

	plotobj1 <- NULL
	plotobj1 <- ggplot(time.weight.sd.data)
	plotobj1 <- plotobj1 + ggtitle("\n(a)\n")
	plotobj1 <- plotobj1 + geom_line(aes(x = TIME,y = SD,colour = Method))
	plotobj1 <- plotobj1 + scale_x_continuous("\nTime since first dose (days)",breaks = c(0,98,210,378))
	plotobj1 <- plotobj1 + scale_y_continuous("Standard deviation of posterior distribution\n")
	plotobj1 <- plotobj1 + theme(legend.position = "none")
	print(plotobj1,vp = vplayout(1,1))

	NTimeWeight <- 1^(TIME.seq-378)
	Peck1005 <- 1.005^(TIME.seq-378)
	Peck101 <- 1.01^(TIME.seq-378)

	time.weight.data <- data.frame(TIME = TIME.seq,
		RELCONT = c(NTimeWeight,Peck1005,Peck101),
		Method = c(rep("NTimeWeight",times = length(TIME.seq)),rep("Peck1.005",times = length(TIME.seq)),rep("Peck1.01",times = length(TIME.seq)))
	)

	plotobj2 <- NULL
	plotobj2 <- ggplot(time.weight.data)
	plotobj2 <- plotobj2 + ggtitle("\n(b)\n")
	plotobj2 <- plotobj2 + geom_line(aes(x = TIME,y = RELCONT,colour = Method))
	plotobj2 <- plotobj2 + scale_x_continuous("\nTime since first dose (days)",breaks = c(0,98,210,378))
	plotobj2 <- plotobj2 + scale_y_continuous("Relative weight in objective function value\n",breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5),labels = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5))
	plotobj2 <- plotobj2 + theme(legend.position = "none")
	print(plotobj2,vp = vplayout(1,2))

	n <- 1000000	# Number of numbers generated from distribution

	# Time == 378
		data1 <- time.weight.sd.data[time.weight.sd.data$TIME == 378,]
		data1.NTimeWeight <- rnorm(n,mean = 0,sd = data1$SD[data1$Method == "NTimeWeight"])
		data1.Peck1005 <- rnorm(n,mean = 0,sd = data1$SD[data1$Method == "Peck1.005"])
		data1.Peck101 <- rnorm(n,mean = 0,sd = data1$SD[data1$Method == "Peck1.01"])

		data1.dist <- data.frame(Method = c(rep("NTimeWeight",times = n),rep("Peck1.005",times = n),rep("Peck1.01",times = n)),
			dist = c(data1.NTimeWeight,data1.Peck1005,data1.Peck101)
		)

		plotobj3 <- NULL
		plotobj3 <- ggplot(data1.dist)
		plotobj3 <- plotobj3 + ggtitle("\n(c)\n")
		plotobj3 <- plotobj3 + geom_density(aes(x = dist,y = ..density..,colour = Method))
		plotobj3 <- plotobj3 + scale_x_continuous("\nPosterior distribution for most recent sample",lim = c(-2,2))
		plotobj3 <- plotobj3 + scale_y_continuous("Density\n",lim = c(0,NA),breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5),labels = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5))
		plotobj3 <- plotobj3 + theme(legend.position = "none")
		print(plotobj3,vp = vplayout(2,1))

	# Time == 210
		data2 <- time.weight.sd.data[time.weight.sd.data$TIME == 210,]
		data2.NTimeWeight <- rnorm(n,mean = 0,sd = data2$SD[data2$Method == "NTimeWeight"])
		data2.Peck1005 <- rnorm(n,mean = 0,sd = data2$SD[data2$Method == "Peck1.005"])
		data2.Peck101 <- rnorm(n,mean = 0,sd = data2$SD[data2$Method == "Peck1.01"])

		data2.dist <- data.frame(Method = c(rep("NTimeWeight",times = n),rep("Peck1.005",times = n),rep("Peck1.01",times = n)),
			dist = c(data2.NTimeWeight,data2.Peck1005,data2.Peck101)
		)

		plotobj4 <- NULL
		plotobj4 <- ggplot(data2.dist)
		plotobj4 <- plotobj4 + ggtitle("\n(d)\n")
		plotobj4 <- plotobj4 + geom_density(aes(x = dist,y = ..density..,colour = Method))
		plotobj4 <- plotobj4 + scale_x_continuous("\nPosterior distribution for sample from 168 days ago",lim = c(-10,10))
		plotobj4 <- plotobj4 + scale_y_continuous("Density\n",lim = c(0,NA),breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5),labels = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5))
		plotobj4 <- plotobj4 + theme(legend.position = "none")
		print(plotobj4,vp = vplayout(2,2))

	# Time == 98
		data3 <- time.weight.sd.data[time.weight.sd.data$TIME == 98,]
		data3.NTimeWeight <- rnorm(n,mean = 0,sd = data3$SD[data3$Method == "NTimeWeight"])
		data3.Peck1005 <- rnorm(n,mean = 0,sd = data3$SD[data3$Method == "Peck1.005"])
		data3.Peck101 <- rnorm(n,mean = 0,sd = data3$SD[data3$Method == "Peck1.01"])

		data3.dist <- data.frame(Method = c(rep("NTimeWeight",times = n),rep("Peck1.005",times = n),rep("Peck1.01",times = n)),
			dist = c(data3.NTimeWeight,data3.Peck1005,data3.Peck101)
		)

		plotobj5 <- NULL
		plotobj5 <- ggplot(data3.dist)
		plotobj5 <- plotobj5 + ggtitle("\n(e)\n")
		plotobj5 <- plotobj5 + geom_density(aes(x = dist,y = ..density..,colour = Method))
		plotobj5 <- plotobj5 + scale_x_continuous("\nPosterior distribution for sample from 280 days ago",lim = c(-20,20))
		plotobj5 <- plotobj5 + scale_y_continuous("Density\n",lim = c(0,NA),breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5),labels = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5))
		plotobj5 <- plotobj5 + theme(legend.position = "none")
		print(plotobj5,vp = vplayout(3,1))

		dev.off()

# ------------------------------------------------------------------------------
# Individual's covariate values changing over time
# Read in population data set
	project.dir <- "/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-bayes/Project/"
	n <- 20  # Number of seed individuals that were simulated
	nsim <- 500 # Number of simulations of seed individuals
	error.model <- "_fourier-extreme-albumin"  # Residual error model (testing proportional and log-normal)
	sim.output.dir <- paste0("SIM",nsim,"_IND",n,error.model)
	output.dir <- paste0("/Volumes/Prosecutor/PhD/InfliximabBayes/infliximab-output/",sim.output.dir,"/")
	pop.data <- read.csv(file = paste0(output.dir,"population_characteristics.csv"))

# Use SIM == 100
	sim.data <- pop.data[pop.data$SIM == 0,]

# Calculate PK parameters as a function of population values, covariates and random effects
 	sim.data$CL <- 0.294*((sim.data$ALB/4)^-1.17)*(1+sim.data$ADA*0.257)*exp(sim.data$ETA1)
	sim.data$V1 <- 3.33*exp(sim.data$ETA2)
	sim.data$Q <- 0.0719*exp(sim.data$ETA3)
	sim.data$V2 <- 1.14*exp(sim.data$ETA4)

# Plot ALB, ADA and ETAs over time
	plotobj6 <- NULL
	plotobj6 <- ggplot(sim.data)
	plotobj6 <- plotobj6 + geom_line(aes(x = TIME,y = ALB,colour = "Albumin (U/L) "))
	plotobj6 <- plotobj6 + geom_step(aes(x = TIME,y = ADA,colour = "ADA Status (0 or 1)  "))
	plotobj6 <- plotobj6 + geom_line(aes(x = TIME,y = CL,colour = "CL (L/day/70 kg) "))
	plotobj6 <- plotobj6 + geom_line(aes(x = TIME,y = V1,colour = "V1 (L/70 kg) "))
	plotobj6 <- plotobj6 + geom_line(aes(x = TIME,y = Q,colour = "Q (L/day/70 kg) "))
	plotobj6 <- plotobj6 + geom_line(aes(x = TIME,y = V2,colour = "V2 (L/70 kg) "))
	plotobj6 <- plotobj6 + scale_x_continuous("\nTime (days)",breaks = seq(0,546,182))
	plotobj6 <- plotobj6 + scale_y_continuous("Variable\n",breaks = c(0,2,4,6,8),labels = c(0,2,4,6,8))
	plotobj6 <- plotobj6 + theme(legend.title = element_blank(),legend.position = "bottom")
	plotobj6 <- plotobj6 + facet_wrap(~ID,ncol = 5,scales = "free_y")
	print(plotobj6)

	ggsave(plot = plotobj6,file = paste0(plotoutput.dir,"seed_characteristics.png"),width = 30,height = 20,units = "cm")
