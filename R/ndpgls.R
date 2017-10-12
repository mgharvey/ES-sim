ndpgls <- function(phy, trait, nd) {
	
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
	colnames(dframe) <- c("species", "trait", "nd")
	data <- comparative.data(data=dframe, phy=phy, names.col="species")

	# PGLS of correlation between inverse splits statistic and trait using Caper

	res <- pgls(nd ~ trait, data=data)

	corr <- summary(res)$coefficients[2,1]
	pval <- summary(res)$coefficients[2,4]

	result <- as.vector(c(corr, pval))
	names(result) <- c("PGLS Slope", "P Value")
	return(result)

}
