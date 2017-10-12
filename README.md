# ES-sim
Code used for testing trait-dependent diversification using tip rate correlation (TRC). 


INTRODUCTION
-------

Code used for testing trait-dependent diversification using tip rate correlation (TRC). 

LICENSE
-------

The code within this repository is available under a 3-clause BSD license. See the License.txt file for more information.

CITATION
--------

If you use this pipeline for your own research, please cite:

* Harvey, MG and Rabosky, DL. In revision. **Continuous traits and speciation rates: Alternatives sto state-dependent diversification models**. 
    
You can also provide a link to this repository if desired:

    https://github.com/mgharvey/ES-sim

USAGE
--------

Requirements:
ape (https://cran.r-project.org/web/packages/ape/)
mvtnorm (https://cran.r-project.org/web/packages/mvtnorm/)

For now, simply download essim.R, then load it in R:

```
source("essim.R")
```

Run it using the command:

```
essim(phy, trait, nsim = 1000)
```

where "phy" is your phylogeny, "trait" is your trait information (with names that match the names on the tips of your phylogeny, and "nsim" is the number of simulations used to build the null distribution of trait-speciation associations for significance testing.

Additional scripts used for simulating and analyzing the data in Harvey and Rabosky (in revision) are included in the directory ./R/ and the simulated datasets examined are in the directory ./data/.