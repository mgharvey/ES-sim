ndsim <- function(phy, trait, nsim = 1000, nd) {
	
	require(ape)
	require(caper)
	require(phytools)
	require(geiger)

	if(missing(nd)) { # If node depth statistics not provided
		
		# Calculate node depth statistic
		rootnode <- length(phy$tip.label) + 1
		nd <- numeric(length(phy$tip.label)) # Number of tips
		for (i in 1:length(nd)){ # For each tip
			node <- i
			nodecount <- 0
			while (node != rootnode){
				node <- phy$edge[,1][phy$edge[,2] == node] # Work back a node	
				nodecount <- nodecount + 1
			}
			nd[i] <- nodecount/max(branching.times(phy))
		}		
	names(nd) <- phy$tip.label
			
	}
		
	# Make phylo comparative data object with trait and inverse splits stat for each species
	dframe <- data.frame(names(trait), trait, log(nd[as.vector(names(trait))]))
	colnames(dframe) <- c("species", "trait", "invsplits")
	data <- comparative.data(data=dframe, phy=phy, names.col="species")

	# Calculate q (rate) using Geiger

	q.trait <- fitContinuous(data$phy, trait, model="BM")
	rate <- q.trait$opt$sigsq
	root <- q.trait$opt$z0

	# Correlation test between inverse splits statistic and trait

	res <- cor.test(data$data$invsplits, data$data$trait, method="pearson")

	# Simulations using Phytools

	sims <- fastBM(data$phy, a = root, sig2 = rate, nsim=nsim) 
	# Root state shouldn't matter, q from empirical data

	# correlation of simulated data

	simdf <- data.frame(rownames(sims), log(nd[as.vector(rownames(sims))]), sims[,1:nsim])
	# Make sure above are all in same order
	colnames(simdf)[1] <- "species"
	colnames(simdf)[2] <- "invsplits"
	sim.data <- comparative.data(data=simdf, phy=phy, names.col="species")
	vars <- sim.data$data[,2:(nsim+1)]
	sim.res <- sapply(vars, function(x) cor.test(sim.data$data$invsplits, x, method="pearson"))
	sim.r <- sim.res[4,]
	sim.p <- sim.res[3,]

	corr <- res$estimate
	upper <- length(sim.r[sim.r >= corr])/nsim
	lower <- length(sim.r[sim.r <= corr])/nsim
	pval <- 2*min(c(upper,lower)) # Two-tailed pval, remove "2" for one-tailed

	result <- as.vector(c(corr, pval))
	names(result) <- c("rho", "P Value")
	return(result)

}
