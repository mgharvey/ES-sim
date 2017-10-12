essim_kendall <- function(phy, trait, nsim = 1000, is) {
	
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
				qx <- qx + el* (1 / 2^(index-1))			
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
	
	# Pearson's correlation between splits statistic and trait
	res <- cor.test(data$data$invsplits, data$data$trait, method="kendall")

	# Brownian simulations 
	vv <- vcv.phylo(as.phylo(data$phy))
	sims <- t(rmvnorm(nsim, sigma=rate*vv))
	rownames(sims) <- rownames(vv)
		
	# Kendall's correlations of simulated datasets
	sim.r <- sapply(1:nsim, function(x) cor.test(log(is[as.vector(rownames(sims))]), sims[,x], method="kendall")$estimate)
	
	# Calculate the two-tailed p value
	corr <- res$estimate
	upper <- length(sim.r[sim.r >= corr])/nsim
	lower <- length(sim.r[sim.r <= corr])/nsim
	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed

	result <- as.vector(c(corr, pval))
	names(result) <- c("rho", "P Value")
	return(result)

}
?sapply