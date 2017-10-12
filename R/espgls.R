espgls <- function(phy, trait, a = 0.5, is) {
	
	require(ape)
	require(caper)
	require(phytools)
	require(geiger)

	if(missing(is)) { # If inverse splits statistics not provided
		
		# Calculate inverse splits statistic
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

	# PGLS of correlation between inverse splits statistic and trait using Caper

	res <- pgls(invsplits ~ trait, data=data)

	corr <- summary(res)$coefficients[2,1]
	pval <- summary(res)$coefficients[2,4]

	result <- as.vector(c(corr, pval))
	names(result) <- c("PGLS Slope", "P Value")
	return(result)

}
