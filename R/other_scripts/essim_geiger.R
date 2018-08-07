essim_geiger <- function(phy, trait, nsim = 1000, es, return.es=FALSE) {
	
	require(ape)
	require(mvtnorm)
	require(caper)
	require(geiger)
	
	#phy <- subtree
	#trait <- trait
	#nsim <- 1000
	#es <- dx
	#return.es=FALSE
	
	if(missing(es)) { # If inverse equal splits statistics not provided, calculate it
		rootnode <- length(phy$tip.label) + 1
		es <- numeric(length(phy$tip.label))
		for (i in 1:length(es)){
			node <- i
			index <- 1
			qx <- 0
			while (node != rootnode){
				el <- phy$edge.length[phy$edge[,2] == node]
				node <- phy$edge[,1][phy$edge[,2] == node]			
				qx <- qx + el* (1 / 2^(index-1))			
				index <- index + 1
			}
			es[i] <- 1/qx
		}		
	names(es) <- phy$tip.label
	}
	
	es <- log(es[phy$tip.label]) # log transform
	trait <- trait[phy$tip.label]
	
	# Pearson's correlation between log inverse equal splits statistic and trait
	res <- cor.test(es, trait, method="pearson")

	# Make phylo comparative data object with trait and inverse splits stat for each species
	dframe <- data.frame(names(trait), trait, es[as.vector(names(trait))])
	colnames(dframe) <- c("species", "trait", "invsplits")
	data <- comparative.data(data=dframe, phy=phy, names.col="species")

	# Fit Brownian motion model to get diffusion rate and root state estimates using GEIGER
	q.trait <- fitContinuous(data$phy, trait, model="BM")
	rate <- q.trait$opt$sigsq
	root <- q.trait$opt$z0
	
	# Brownian simulations 
	sims <- t(rmvnorm(nsim, sigma=rate*vv))
	rownames(sims) <- rownames(vv)
		
	# Pearson's correlations of simulated datasets
	sim.r <- sapply(1:nsim, function(x) cor.test(es[as.vector(rownames(sims))], sims[,x], method="pearson")$estimate)
	
	# Calculate the two-tailed p value
	corr <- res$estimate
	upper <- (length(sim.r[sim.r >= corr])+1)/(nsim+1)
	lower <- (length(sim.r[sim.r <= corr])+1)/(nsim+1)
	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed

	if(missing(return.es)) { # output just rho and p value
		result <- as.vector(c(corr, pval))
		names(result) <- c("rho", "P Value")
		return(result)
	} else { # output rho, p value, and list of es values
		result <- as.vector(c(corr, pval, list(es)))
		names(result) <- c("rho", "P Value", "es")
		return(result)		
	}

}
