ndsim <- function(phy, trait, nsim = 1000, nd) {
	
	require(ape)
	require(mvtnorm)

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
		
	nd <- log(nd[phy$tip.label]) # log transform
	trait <- trait[phy$tip.label]
	
	# Pearson's correlation between splits statistic and trait
	res <- cor.test(nd, trait, method="pearson")

	# Fit Brownian motion model to get diffusion rate and root state estimates
	vv <- vcv.phylo(as.phylo(phy))
	onev <- matrix(rep(1, length(trait)), nrow=length(trait), ncol=1)
	root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% trait))
	rate <- as.vector((t(trait-root) %*% solve(vv) %*% (trait-root))/length(trait))
	
	# Brownian simulations 
	sims <- t(rmvnorm(nsim, sigma=rate*vv))
	rownames(sims) <- rownames(vv)
		
	# Pearson's correlations of simulated datasets
	sim.r <- sapply(1:nsim, function(x) cor.test(nd[as.vector(rownames(sims))], sims[,x], method="pearson")$estimate)
	
	# Calculate the two-tailed p value
	corr <- res$estimate
	upper <- length(sim.r[sim.r >= corr])/nsim
	lower <- length(sim.r[sim.r <= corr])/nsim
	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed

	result <- as.vector(c(corr, pval))
	names(result) <- c("rho", "P Value")
	return(result)

}
