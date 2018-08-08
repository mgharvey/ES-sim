require(ape);
require(mvtnorm);

essim <- function(phy, trait, nsim = 1000, es, return.es=FALSE) {
	
    if (missing(es)) { # If inverse equal splits statistics not provided, calculate it
        #rootnode <- length(phy$tip.label) + 1
        #es <- numeric(length(phy$tip.label))
        #for (i in 1:length(es)){
        #	node <- i
        #	index <- 1
        #	qx <- 0
        #	while (node != rootnode){
        #		el <- phy$edge.length[phy$edge[,2] == node]
        #		node <- phy$edge[,1][phy$edge[,2] == node]			
        #		qx <- qx + el* (1 / 2^(index-1))			
        #		index <- index + 1
        #	}
        #	es[i] <- 1/qx
        #}		
        #names(es) <- phy$tip.label
        es <- compute_es(phy);
    }
    
    # check if all tips have traits. will work whether newly computed or passed in
    if (length(trait) < length(phy$tip.label)) {
        idx <- which(!phy$tip.label %in% names(trait));
        phy <- drop.tip(phy, idx);
    }
	
	es <- log(es[phy$tip.label]) # log transform
	trait <- trait[phy$tip.label]
	
	# Pearson's correlation between log inverse equal splits statistic and trait
	res <- cor.test(es, trait, method="pearson")

	# Fit Brownian motion model to get diffusion rate and root state estimates
	vv <- vcv.phylo(as.phylo(phy))
	onev <- matrix(rep(1, length(trait)), nrow=length(trait), ncol=1)
	root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% trait))
	rate <- as.vector((t(trait-root) %*% solve(vv) %*% (trait-root))/length(trait))
	
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

require(phangorn); # for Ancestors
compute_es <- function (phy) {
    Ntip <- length(phy$tip.label);
    rootnd <- Ntip + 1L;
    
    es <- numeric(Ntip);
    
    for (k in 1:Ntip) {
        # get lineage
        lin <- c(k, Ancestors(phy, k));
        # drop root node
        lin <- lin[-length(lin)];
        
        # set indices
        inds <- seq(1, length(lin), 1);
        # get els
        els <- phy$edge.length[match(lin, phy$edge[,2])];
        # compute in one
        es[k] <- 1/sum(els * (1 / 2^(inds-1)));
        
    }
    names(es) <- phy$tip.label;
    return (es);
}
