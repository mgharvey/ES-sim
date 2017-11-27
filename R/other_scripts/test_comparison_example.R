library(ape)
library(mvtnorm)
library(diversitree)
source("./R/essim.R")

# Power test

trees250sdd <- read.tree("./data/trees_Core/trees250sdd.txt")
states250sdd <- read.table("./data/traits_Core/trees250sdd_states.txt", header=TRUE, row.names=1)

essim.p <- vector()
q.p <- vector()

for (i in 1:length(trees250sdd)) {

	print(i)
	trait <- states250sdd[,i]
	names(trait) <- row.names(states250sdd)
	
	# ES-sim

	essim.res <- essim(trees250sdd[[i]], trait, nsim = 1000)
	essim.p <- c(essim.p, essim.res[2])
		
	# QuaSSE

	p.constant <- starting.point.quasse(trees250sdd[[i]], trait)
	xr <- range(trait) + c(-1,1) * 20 * p.constant["diffusion"]
	linear.x <- make.linear.x(xr[1], xr[2])
	control <- list(parscale=.1, reltol=0.001)
	nodrift <- function(f) constrain(f, drift ~ 0)
	lik.constant <- make.quasse(trees250sdd[[i]], trait, (1/200), constant.x, constant.x)
	fit.constant <- find.mle(nodrift(lik.constant), p.constant, lower=0, control=control, verbose=0)			 
	p <- c(fit.constant$par[1], l.m.=0, fit.constant$par[2:3])
	lik <- make.quasse(trees250sdd[[i]], trait, (1/200), linear.x, constant.x)
	fit <- find.mle(nodrift(lik), p, control=control, verbose=0)			
	sim.res.q <- anova(fit.constant, fit)
	q.p <- c(q.p, sim.res.q[2,5])

}

essim.pow <- length(essim.p[essim.p < 0.05])/length(trees250sdd)
quasse.pow <- length(q.p[q.p < 0.05])/length(trees250sdd)


# FDR test

trees250nsdd <- read.tree("./data/trees_Core/trees250nsdd.txt")
states250nsdd <- read.table("./data/traits_Core/trees250nsdd_states.txt", header=TRUE, row.names=1)

essim.p <- vector()
q.p <- vector()

for (i in 1:length(trees250nsdd)) {
	
	print(i)
	trait <- states250nsdd[,i]
	names(trait) <- row.names(states250nsdd)
	
	# ES-sim

	essim.res <- essim(trees250nsdd[[i]], trait, nsim = 1000)
	essim.p <- c(essim.p, essim.res[2])
		
	# QuaSSE

	p.constant <- starting.point.quasse(trees250nsdd[[i]], trait)
	xr <- range(trait) + c(-1,1) * 20 * p.constant["diffusion"]
	linear.x <- make.linear.x(xr[1], xr[2])
	control <- list(parscale=.1, reltol=0.001)
	nodrift <- function(f) constrain(f, drift ~ 0)
	lik.constant <- make.quasse(trees250nsdd[[i]], trait, (1/200), constant.x, constant.x)
	fit.constant <- find.mle(nodrift(lik.constant), p.constant, lower=0, control=control, verbose=0)			 
	p <- c(fit.constant$par[1], l.m.=0, fit.constant$par[2:3])
	lik <- make.quasse(trees250nsdd[[i]], trait, (1/200), linear.x, constant.x)
	fit <- find.mle(nodrift(lik), p, control=control, verbose=0)			
	sim.res.q <- anova(fit.constant, fit)
	q.p <- c(q.p, sim.res.q[2,5])

}

essim.fdr <- length(essim.p[essim.p < 0.05])/length(trees250nsdd)
quasse.fdr <- length(q.p[q.p < 0.05])/length(trees250nsdd)


# FDR with background diversification dynamics test

for (i in 1:length(trees250sdd)) {

	print(i)
	oldtrait <- states250sdd[,i]
	names(oldtrait) <- row.names(states250sdd)

	# Simulate a new trait with brownian motion

	# Make phylo comparative data object with trait and inverse splits stat for each species
	dframe <- data.frame(names(oldtrait), oldtrait)
	colnames(dframe) <- c("species", "oldtrait")
	data <- comparative.data(data=dframe, phy=trees250sdd[[i]], names.col="species")
	q.trait <- fitContinuous(data$phy, oldtrait, model="BM")
	rate <- q.trait$opt$sigsq
	root <- q.trait$opt$z0
	trait <- sim.character(trees250sdd[[i]], rate, x0=root, model="bm")
	
	# ES-sim

	essim.res <- essim(trees250sdd[[i]], trait, nsim = 1000)
	essim.p <- c(essim.p, essim.res[2])
		
	# QuaSSE

	p.constant <- starting.point.quasse(trees250sdd[[i]], trait)
	xr <- range(trait) + c(-1,1) * 20 * p.constant["diffusion"]
	linear.x <- make.linear.x(xr[1], xr[2])
	control <- list(parscale=.1, reltol=0.001)
	nodrift <- function(f) constrain(f, drift ~ 0)
	lik.constant <- make.quasse(trees250sdd[[i]], trait, (1/200), constant.x, constant.x)
	fit.constant <- find.mle(nodrift(lik.constant), p.constant, lower=0, control=control, verbose=0)			 
	p <- c(fit.constant$par[1], l.m.=0, fit.constant$par[2:3])
	lik <- make.quasse(trees250sdd[[i]], trait, (1/200), linear.x, constant.x)
	fit <- find.mle(nodrift(lik), p, control=control, verbose=0)			
	sim.res.q <- anova(fit.constant, fit)
	q.p <- c(q.p, sim.res.q[2,5])	

}

essim.fdrb <- length(essim.p[essim.p < 0.05])/length(trees250nsdd)
quasse.fdrb <- length(q.p[q.p < 0.05])/length(trees250nsdd)

