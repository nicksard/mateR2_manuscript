# Not sourced. Exists so renv's dependency scan sees packages that are used
# at runtime through mateR2's Suggests rather than called directly here.
library(coda)    # run_mcmc_chains() returns coda::mcmc.list traces
library(igraph)  # count_network_components() / comating_pair_density()
