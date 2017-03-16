# OncoSimulR_NKmodel
OncoSimulR package w. NKmodel fitness
Package: OncoSimulR
Type: Package
Title: Forward Genetic Simulation of Cancer Progression with Epistasis 
Version: 2.5.12
Date: 2017-02-08
Authors@R: c(person("Ramon", "Diaz-Uriarte", role = c("aut", "cre"),
		     email = "rdiaz02@gmail.com"),
	      person("Mark", "Taylor", role = "ctb", email = "ningkiling@gmail.com"))
Author: Ramon Diaz-Uriarte [aut, cre],
	Mark Taylor [ctb]
Maintainer: Ramon Diaz-Uriarte <rdiaz02@gmail.com>
Description: Functions for forward population genetic simulation in
    asexual populations, with special focus on cancer progression.
    Fitness can be an arbitrary function of genetic interactions between
    multiple genes or modules of genes, including epistasis, order
    restrictions in mutation accumulation, and order effects.  Mutation
    rates can differ between genes, and we can include mutator/antimutator
    genes (to model mutator phenotypes). Simulations
    use continuous-time models and can include driver and passenger genes
    and modules.  Also included are functions for: simulating random DAGs
    of the type found in Oncogenetic Tress, Conjunctive Bayesian Networks,
    and other tumor progression models; plotting and sampling from
    single or multiple realizations of the simulations, including
    single-cell sampling; plotting the parent-child relationships of the
    clones; generating random fitness landscapes (Rough Mount Fuji, House
    of Cards, and additive models) and plotting them.
biocViews: BiologicalQuestion, SomaticMutation
License: GPL (>= 3)
URL: https://github.com/rdiaz02/OncoSimul, https://popmodels.cancercontrol.cancer.gov/gsr/packages/oncosimulr/
BugReports: https://github.com/rdiaz02/OncoSimul/issues
Depends: R (>= 3.3.0)
Imports: Rcpp (>= 0.12.4), parallel, data.table, graph, Rgraphviz, gtools, igraph, methods, RColorBrewer, grDevices, car, dplyr, smatr, ggplot2, ggrepel, nem
Suggests: BiocStyle, knitr, Oncotree, testthat (>= 1.0.0), rmarkdown, bookdown, pander
LinkingTo: Rcpp
VignetteBuilder: knitr
RoxygenNote: 6.0.1.9000
