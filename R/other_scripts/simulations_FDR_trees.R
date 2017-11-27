setwd("/Users/michaelharvey/Documents/research/tdd/fdr_scenarios/")
getwd()

# Script: 
# By: Michael G. Harvey 
# Date:

# A script to simulate trees.

library(ape)
library(geiger)
library(TreeSim)
library(diversitree)

### SIMULATE TREES ###

# Constant rate 

trees.constant <- vector("list", 250) # Vector to place trees
class(trees.constant) <- "multiPhylo"
t = 0
while(t < 250) {
	ntax <- 200 # number of tips in tree
	tree <- sim.bd.age(3, 1, 2, 0.2, complete=FALSE)[[1]]
	try(if(length(tree) == 5) { # Remove null trees produced by TreeSim
		if(length(tree$tip.label) < 250) {
			if(length(tree$tip.label) >= 100) {
				t = t+1		
				trees.constant[[t]] <- tree
				print(t)
			}
		}
	})
}
write.tree(trees.constant, "trees_constant.txt")

# Slowdown

trees.slowdown <- vector("list", 250) # Vector to place trees
class(trees.slowdown) <- "multiPhylo"
t = 0
while(t < 250) {
	tree <- sim.bd.age(5, 1, 2, 0.2, complete=FALSE, K=300)[[1]]
	try(if(length(tree) == 5) { # Remove null trees produced by TreeSim
		if(length(tree$tip.label) < 250) {
			if(length(tree$tip.label) >= 100) {
				t = t+1		
				trees.slowdown[[t]] <- tree	
				print(t)	
			}
		}		
	})
}
write.tree(trees.slowdown, "trees_slowdown.txt")

# True QuaSSE

trees.quasse <- vector("list", 250) # Vector to place trees
class(trees.quasse) <- "multiPhylo"
t = 0
while(t < 250) {
	m <- 0.004 # strength of relationship between trait and speciation
	ntax <- 200 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 0.06)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=200, x0=0, single.lineage=FALSE, verbose=FALSE)
	t = t+1
	trees.quasse[[t]] <- phy
	print(t)
}
write.tree(trees.quasse, "trees_quasse.txt")

# True BiSSE

trees.bisse <- vector("list", 250) # Vector to place trees
class(trees.quasse) <- "multiPhylo"
t = 0
while(t < 250) {
	m <- 0.004 # strength of relationship between trait and speciation
	ntax <- 200 # number of tips in tree
	phy <- tree.bisse(c(0.1, 0.3, 0.03, 0.03, 0.01, 0.01), max.taxa=200, x0=0)
	t = t+1
	trees.bisse[[t]] <- phy
	print(t)
}
write.tree(trees.quasse, "trees_bisse.txt")

# Empirical 1 (coral supertree) - taken from Dryad tree set AJW (Rabosky and Goldberg 2017)

# Empirical 2 (carnivore tree) - taken from Dryad tree QVW (Rabosky and Goldberg 2017), only 1 tree

# BAMM trees ddK4 - taken from Dryad (Rabosky 2014)





