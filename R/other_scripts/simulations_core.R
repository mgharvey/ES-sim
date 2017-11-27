library(diversitree)

### Representative set of 100 simulations with and without trait-dependence ###

# 50 taxa, SDD (state-dependent diversification)

trees50sdd <- vector("list", 100) # Vector to place trees
class(trees50sdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {	
	m <- 0.004 # strength of relationship between trait and speciation
	ntax <- 50 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 0.06)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees50sdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees50sdd, "trees50sdd.txt")
write.table(as.matrix(states), "trees50sdd_states.txt")

# 50 taxa, non-SDD

trees50nsdd <- vector("list", 100) # Vector to place trees
class(trees50nsdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0 # strength of relationship between trait and speciation
	ntax <- 50 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 0.06) 
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees50nsdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees50nsdd, "trees50nsdd.txt")
write.table(as.matrix(states), "trees50nsdd_states.txt")

### The full set of simulations examined, tweaking different parameters ###

# 250 taxa, SDD

trees250sdd <- vector("list", 100) # Vector to place trees
class(trees250sdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0.004 # strength of relationship between trait and speciation
	ntax <- 250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 0.06) 
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees250sdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees250sdd, "trees250sdd.txt")
write.table(as.matrix(states), "trees250sdd_states.txt")

# 250 taxa, non-SDD

trees250nsdd <- vector("list", 100) # Vector to place trees
class(trees250nsdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0 # strength of relationship between trait and speciation
	ntax <- 250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 0.06)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees250nsdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees250nsdd, "trees250nsdd.txt")
write.table(as.matrix(states), "trees250nsdd_states.txt")

# 1250 taxa, SDD

trees1250sdd <- vector("list", 100) # Vector to place trees
class(trees1250sdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {	
	m <- 0.004 # strength of relationship between trait and speciation
	ntax <- 1250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 0.06)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees1250sdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees1250sdd, "trees1250sdd.txt")
write.table(as.matrix(states), "trees1250sdd_states.txt")

# 1250 taxa, non-SDD

trees1250nsdd <- vector("list", 100) # Vector to place trees
class(trees1250nsdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0 # strength of relationship between trait and speciation
	ntax <- 1250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 0.06)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees1250nsdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees1250nsdd, "trees1250nsdd.txt")
write.table(as.matrix(states), "trees1250nsdd_states.txt")

# slow trait change

trees250sdd <- vector("list", 100) # Vector to place trees
class(trees250sdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0.004 # strength of relationship between trait and speciation
	ntax <- 250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 0.006) 
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees250sdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees250sdd, "trees250sddslow.txt")
write.table(as.matrix(states), "trees250sddslow_states.txt")

# very slow trait change

trees250nsdd <- vector("list", 100) # Vector to place trees
class(trees250nsdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0 # strength of relationship between trait and speciation
	ntax <- 250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 0.0006)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees250nsdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees250nsdd, "trees250nsddveryslow.txt")
write.table(as.matrix(states), "trees250nsddveryslow_states.txt")

# very, very slow trait change

trees250nsdd <- vector("list", 100) # Vector to place trees
class(trees250nsdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0 # strength of relationship between trait and speciation
	ntax <- 250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 0.00006)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees250nsdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees250nsdd, "trees250nsddveryveryslow.txt")
write.table(as.matrix(states), "trees250nsddveryveryslow_states.txt")

# fast trait change

trees250nsdd <- vector("list", 100) # Vector to place trees
class(trees250nsdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0 # strength of relationship between trait and speciation
	ntax <- 250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 0.6)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees250nsdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees250nsdd, "trees250nsddfast.txt")
write.table(as.matrix(states), "trees250nsddfast_states.txt")

# very fast trait change

trees250nsdd <- vector("list", 100) # Vector to place trees
class(trees250nsdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0 # strength of relationship between trait and speciation
	ntax <- 250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 6)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees250nsdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees250nsdd, "trees250nsddveryfast.txt")
write.table(as.matrix(states), "trees250nsddveryfast_states.txt")

# very, very fast trait change

trees250nsdd <- vector("list", 100) # Vector to place trees
class(trees250nsdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0 # strength of relationship between trait and speciation
	ntax <- 250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	char <- make.brownian.with.drift(0, 60)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees250nsdd[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees250nsdd, "trees250nsddveryveryfast.txt")
write.table(as.matrix(states), "trees250nsddveryveryfast_states.txt")

# OU trait evolution, very weak "pull"

trees <- vector("list", 100) # Vector to place trees
class(trees) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0.004 # strength of relationship between trait and speciation
	ntax <- 250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	make.ou.mgh <- function (drift, diffusion, a, theta) 
	function(x, dt) x + rnorm(length(x), (a*(theta-x)) * dt, sqrt(dt * diffusion))
	char <- make.ou.mgh(0, 0.06, 0.002, 2)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees, "trees250sdd_OUweak.txt")
write.table(as.matrix(states), "trees250sdd_OUweak_states.txt")

# OU trait evolution, weak "pull"

trees <- vector("list", 100) # Vector to place trees
class(trees) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0.004 # strength of relationship between trait and speciation
	ntax <- 250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	make.ou.mgh <- function (drift, diffusion, a, theta) 
	function(x, dt) x + rnorm(length(x), (a*(theta-x)) * dt, sqrt(dt * diffusion))
	char <- make.ou.mgh(0, 0.06, 0.02, 2)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees, "trees250sdd_OUweak.txt")
write.table(as.matrix(states), "trees250sdd_OUweak_states.txt")

# OU trait evolution, strong "pull"

trees <- vector("list", 100) # Vector to place trees
class(trees) <- "multiPhylo"
states <- vector()
for (i in 1:100) {
	m <- 0.004 # strength of relationship between trait and speciation
	ntax <- 250 # number of tips in tree
	linear.x <- make.linear.x(-5, 5)
	lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
	mu <- function(x) constant.x(x, 0)
	make.ou.mgh <- function (drift, diffusion, a, theta) 
	function(x, dt) x + rnorm(length(x), (a*(theta-x)) * dt, sqrt(dt * diffusion))
	char <- make.ou.mgh(0, 0.06, 0.2, 2)
	phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=TRUE)
	trees[[i]] <- phy
	states <- cbind(states, phy$tip.state)
}
write.tree(trees, "trees250sdd_OUweak.txt")
write.table(as.matrix(states), "trees250sdd_OUweak_states.txt")





