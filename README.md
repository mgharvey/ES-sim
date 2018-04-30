INTRODUCTION
-------

Code used for testing trait-dependent diversification using tip rate correlation (TRC). 

LICENSE
-------

The code within this repository is available under a 3-clause BSD license. See the License.txt file for more information.

CITATION
--------

If you use this pipeline for your own research, please cite:

* Harvey, MG and Rabosky, DL. In press. **Continuous traits and speciation rates: Alternatives to state-dependent diversification models**. *Methods in Ecology and Evolution* doi: 10.1111/2041-210X.12949
    
You can also provide a link to this repository if desired:

    https://github.com/mgharvey/ES-sim

USAGE
--------

Code used for simulating and analyzing the data in Harvey and Rabosky (in press) are included in the directory ./R/ and the simulated datasets examined are in the directory ./data/.

Requirements:

* ape (https://cran.r-project.org/web/packages/ape/)

* mvtnorm (https://cran.r-project.org/web/packages/mvtnorm/)

For now, simply download essim.R from the ./R/ directory of this repository, then load it in R:

```
source("essim.R")
```

Run it using the command

```
essim(phy, trait, nsim = 1000)
```

where "phy" is your phylogeny, "trait" is a vector containing your trait information (with names that match the names on the tips of your phylogeny), and "nsim" is the number of simulations used to build the null distribution of trait-speciation associations for significance testing. There are also optional arguments "return.is", which if "True" returns a named list of log inverse equal splits values, and "is", which can be used to supply an existing named vector of log inverse equal splits values. The test assumptions are the same as for fitting a Brownian motion model to phylogenetic comparative data - that your data fit a multivariate normal distribution with the trait covariance between tips determined by the amount of time they have shared a common ancestor. The test will return the Pearson's correlation coefficient (rho) and the simulation-based two-tailed p-value.

DOI
--------

[![DOI](https://zenodo.org/badge/106715524.svg)](https://zenodo.org/badge/latestdoi/106715524)


