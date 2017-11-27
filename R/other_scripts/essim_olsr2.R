essim_olsr2 <- function(phy, trait, a = 0.5, nsim = 1000, is) {
	
	require(ape)
	require(mvtnorm)

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
	
	is <- log(is[phy$tip.label]) # log transform
	trait <- trait[phy$tip.label]

	# OLS of correlation between inverse splits statistic and trait using Caper
	res <- lm(is ~ trait)

	# Fit Brownian motion model to get diffusion rate and root state estimates
	vv <- vcv.phylo(as.phylo(phy))
	onev <- matrix(rep(1, length(trait)), nrow=length(trait), ncol=1)
	root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% trait))
	rate <- as.vector((t(trait-root) %*% solve(vv) %*% (trait-root))/length(trait))

	# Brownian simulations 
	sims <- t(rmvnorm(nsim, sigma=rate*vv))
	rownames(sims) <- rownames(vv)
		
	# OLS of simulated datasets
	simdf <- data.frame(rownames(sims), is[as.vector(rownames(sims))], sims[,1:nsim])
	colnames(simdf)[1] <- "species"
	colnames(simdf)[2] <- "invsplits"
	vars <- simdf[,3:(nsim+2)]
	sim.r <- sapply(vars, function(x) summary(lm(invsplits ~ x, data = simdf))$r.squared)
	
	# Calculate the two-tailed p value
	corr <- summary(res)$r.squared
	upper <- length(sim.r[sim.r >= corr])/nsim
	lower <- length(sim.r[sim.r <= corr])/nsim
	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed

	result <- as.vector(c(corr, pval))
	names(result) <- c("slope", "P Value")
	return(result)

}
