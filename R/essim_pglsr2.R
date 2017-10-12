essim_pglsr2 <- function(phy, trait, a = 0.5, nsim = 1000, is) {
	
	require(ape)
	require(caper)
	require(geiger)

	if(missing(is)) { # If inverse splits statistics not provided, calculate it
		rootnode <- length(phy$tip.label) + 1
		is <- numeric(length(phy$tip.label))
		for (i in 1:length(is)){
			node <- i
			index <- 1
			qx <- 0
			while (node != rootnode){
				el <- phy$edge.length[phy$edge[,2] == node]
				node <- phy$edge[,1][phy$edge[,2] == node]			
				qx <- qx + el* (1 / (1/a)^(index-1))			
				index <- index + 1
			}
			is[i] <- 1/qx
		}		
	names(is) <- phy$tip.label
	}
	
	# Make phylo comparative data object with trait and inverse splits stat for each species
	dframe <- data.frame(names(trait), trait, log(is[as.vector(names(trait))]))
	colnames(dframe) <- c("species", "trait", "invsplits")
	data <- comparative.data(data=dframe, phy=phy, names.col="species")

	# Fit Brownian motion model to get diffusion rate and root state estimates using GEIGER
	q.trait <- fitContinuous(data$phy, trait, model="BM")
	rate <- q.trait$opt$sigsq
	root <- q.trait$opt$z0
	
	# PGLS of correlation between inverse splits statistic and trait using Caper
	res <- pgls(invsplits ~ trait, data=data)

	# Brownian simulations 
	vv <- vcv.phylo(as.phylo(data$phy))
	sims <- t(rmvnorm(nsim, sigma=rate*vv))
	rownames(sims) <- rownames(vv)
		
	# PGLS of simulated datasets
	simdf <- data.frame(rownames(sims), log(is[as.vector(rownames(sims))]), sims[,1:nsim])
	colnames(simdf)[1] <- "species"
	colnames(simdf)[2] <- "invsplits"
	sim.data <- comparative.data(data=simdf, phy=phy, names.col="species")
	vars <- sim.data$data[,2:(nsim+1)]
	sim.r <- sapply(vars, function(x) summary(pgls(invsplits ~ x, data = sim.data))$r.squared)
	
	# Calculate the two-tailed p value
	corr <- summary(res)$r.squared
	upper <- length(sim.r[sim.r >= corr])/nsim
	lower <- length(sim.r[sim.r <= corr])/nsim
	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed

	result <- as.vector(c(corr, pval))
	names(result) <- c("slope", "P Value")
	return(result)

}
