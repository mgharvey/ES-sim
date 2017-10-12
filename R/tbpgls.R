tbpgls <- function(phy, trait) {
	
	require(ape)
	require(caper)
	require(phytools)
	require(geiger)

	# Calculate terminal edge lengths
	n <- length(phy$tip.label)
	# from Liam Revell:
	invis <- setNames(phy$edge.length[sapply(1:n, function(x,y) which(y==x), y=phy$edge[,2])], phy$tip.label)
	tb <- 1/invis
	
	# Make phylo comparative data object with trait and inverse splits stat for each species
	dframe <- data.frame(names(trait), trait, log(tb[as.vector(names(trait))]))
	colnames(dframe) <- c("species", "trait", "tb")
	data <- comparative.data(data=dframe, phy=phy, names.col="species")

	# PGLS of correlation between inverse splits statistic and trait using Caper

	res <- pgls(tb ~ trait, data=data)

	corr <- summary(res)$coefficients[2,1]
	pval <- summary(res)$coefficients[2,4]

	result <- as.vector(c(corr, pval))
	names(result) <- c("PGLS Slope", "P Value")
	return(result)

}
