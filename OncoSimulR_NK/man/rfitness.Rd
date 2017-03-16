\name{rfitness}
\alias{rfitness}


\title{Generate random fitness.}

\description{ Generate random fitness landscapes under a House of Cards,
  Rough Mount Fuji, or additive model.  }


\usage{

rfitness(g, c = 0.5, sd = 1, mu = 1, reference = "random", k = NULL, k_model = c(NA, "arithmetic", "geometric"), scale = NULL,
         wt_is_1 = c("subtract", "divide", "force", "no"),
         log = FALSE, min_accessible_genotypes = NULL,
         accessible_th = 0, truncate_at_0 = TRUE)
}




\arguments{

  \item{g}{Number of genes.}

  \item{c}{The decrease in fitness of a genotype per each unit increase
    in Hamming distance from the reference genotype (see \code{reference}).}

  \item{sd}{The standard deviation of the random component (a normal
  distribution of mean \code{mu} and standard deviation \code{sd}).}

\item{mu}{The mean of the random component (a normal distribution of
mean \code{mu} and standard deviation \code{sd}).}


\item{reference}{The reference genotype: for the deterministic, additive
  part, this is the genotype with maximal fitness, and all other
  genotypes decrease their fitness by \code{c} for every unit of Hamming
  distance from this reference. If "random" a genotype will be randomly
  chosen as the reference. If "max" the genotype with all positions
  mutated will be chosen as the reference. If you pass a vector (e.g.,
  \code{reference = c(1, 0, 1, 0)}) that will be the reference genotype.
  If "random2" a genotype will be randomly chosen as the reference. In
  contrast to "random", however, not all genotypes have the same
  probability of being chosen; here, what is equal is the probability
  that the reference genotype has 1, 2, ..., g, mutations (and, once a
  number mutations is chosen, all genotypes with that number of
  mutations have equal probability of being the reference). }

\item{k}{ Number of neighbours, measures the interaction between genes. To investigate the effect of epistatic interactions on the adaptive process, we can use the NK model, which is a sytem used to study interactions between loci with different alleles. In a NK model, each loci (g) contributes to the fitness of the organism via an interaction with k adjacent loci. \code{k = NULL}}

\item{k_model}{Either NA, "arithmetic" or "geometric". The overall fitness of an organism is usually given by the average of the g fitness components, but you can also use the geometric mean for its calculation, this way, one could then introduce lethal mutations by setting one or more elements to zero. \code{k_model = c(NA, "arithmetic", "geometric")}}

\item{scale}{Either NULL (nothing is done) or a two-element vector. If a
  two-element vector, fitness is re-scaled between \code{scale[1]} (the
  minimum) and \code{scale[2]} (the maximum).}

\item{wt_is_1}{If "divide" (the default) the fitness of all genotypes is
  divided by the fitness of the wildtype (after possibly adding a value
  to ensure no negative fitness) so that the wildtype (the genotype with
  no mutations) has fitness 1. This is a case of scaling, and it is
  applied after \code{scale}, so if you specify both
  "wt_is_1 = 'divide'" and use an argument for \code{scale} it is most
  likely that the final fitness will not respect the limits in
  \code{scale}.

  If "subtract" we shift all the fitness values (subtracting fitness of
  the wildtype and adding 1) so that the wildtype ends up with a fitness
  of 1. This is also applied after \code{scale}, so if you specify both
  "wt_is_1 = 'subtract'" and use an argument for \code{scale} it is most
  likely that the final fitness will not respect the limits in
  \code{scale} (though the distorsion might be simpler to see as just a
  shift up or down).
  
  If "force" we simply set the fitness of the wildtype to 1, without any
  divisions. This means that the \code{scale} argument would work (but
  it is up to you to make sure that the range of the scale argument
  includes 1 or you might not get what you want). Note that using this
  option can easily lead to landscapes with no accessible genotypes
  (unless you also use \code{scale}).

  If "none", the fitness of the wildtype is not touched.  }


\item{log}{If TRUE, log-transform fitness.}

\item{min_accessible_genotypes}{If not NULL, the minimum number of
  accessible genotypes in the fitness landscape. A genotype is
  considered accessible if you can reach if from the wildtype by going
  through at least one path where all changes in fitness are larger or
  equal to \code{accessible_th}. The changes in fitness are considered
  at each mutational step, i.e., at each addition of one mutation we
  compute the difference between the genotype with \code{k + 1}
  mutations minus the ancestor genotype with \code{k} mutations. Thus, a
  genotype is considered accessible if there is at least one path where
  fitness increases at each mutational step by at least
  \code{accessible_th}.

  If the condition is not satisfied, we continue generating random
  fitness landscapes with the specified parameters until the condition
  is satisfied.

  (Why check against NULL and not against zero? Because this allows you
  to count accessible genotypes even if you do not want to ensure a
  minimum number of accessible genotypes.)
}

\item{accessible_th}{The threshold for the minimal change in fitness at
  each mutation step (i.e., between successive genotypes) that allows a
  genotype to be regarded as accessible. This only applies if
  \code{min_accessible_genotypes} is larger than 0.  So if you want to
  allow small decreases in fitness in successive steps, use a small
  negative value for \code{accessible_th}.  }

\item{truncate_at_0}{If TRUE (the default) any fitness <= 0 is
  substituted by a small positive constant (1e-9). Why? Because
  MAGELLAN and some plotting routines can have trouble (specially if you
  log) with values <=0. Or we might have trouble if we want to log the
  fitness.}

} 


\details{

  The model used here follows the Rough Mount Fuji model in Szendro et
  al., 2013 or Franke et al., 2011. Fitness is given as

  \deqn{f(i) = -c d(i, reference) + x_i}

  where \eqn{d(i, j)} is the Hamming distance between genotypes \eqn{i}
  and \eqn{j} (the number of positions that differ) and \eqn{x_i} is a
  random variable (in this case, a normal deviate of mean \code{mu}
  and standard deviation \code{sd}).

  Setting \eqn{c = 0} we obtain a House of Cards model. Setting \eqn{sd
    = 0} fitness is given by the distance from the reference and if the
    reference is the genotype with all positions mutated, then we have a
    fully additive model (fitness increases linearly with the number of
    positions mutated).

  For OncoSimulR, we often want the wildtype to have a mean of
  1. Reasonable settings are \code{mu = 1} and \code{wt_is_1 =
  'subtract'} so that we simulate from a distribution centered in 1, and
  we make sure afterwards (via a simple shift) that the wildtype is
  actuall 1. The \code{sd} controls the standard deviation, with the
  usual working and meaning as in a normal distribution, unless \code{c}
  is different from zero. In this case, with \code{c} large, the range
  of the data can be large, specially if \code{g} (the number of genes)
  is large.
  
  On the other hand, to investigate the epistatic interactions on the adaptive process, we can employ the NK model, which is a simple system used to study interactions between loci with different alleles. The model consists of N loci in circular, binary sequences (Kauffman 1993; Kauffman and Levin 1987). Each of the N (g) loci contributes to the fitness of the organism via an interaction with K (k) adjacent loci. For each locus i a lookup-table consisting of uniform random numbers represents the fitness component of a binary sequence of lengthh K + 1. For example, K = 1 (interaction with one other locus) is modeled by creating random numbers for the four possible binary pairs 00, 01, 10, 11 for each of the g loci, that is, the fitness component at one locus is conditional on the allele at the other one locus. The overall fitness of an organism is usually given by the average of the N fitness components, we can also use the geometric mean (motivated by the fact that one could then introduce lethal mutations by setting one or more elements in the lookup -tables to zero.)
  
} 

\value{
  
  An matrix with \code{g + 1} columns. Each column corresponds to a
  gene, except the last one that corresponds to fitness. 1/0 in a gene
  column denotes gene mutated/not-mutated. (For ease of use in other
  functions, this matrix has class  "genotype_fitness_matrix".) 

  If you have specified \code{min_accessible_genotypes > 0}, the return
  object has added attributes \code{accessible_genotypes} and
  \code{accessible_th} that show the number of accessible
  genotypes under the specified  threshold.
}
  
\references{

  Szendro I.~G. et al. (2013). Quantitative analyses of empirical
fitness landscapes. \emph{Journal of Statistical Mehcanics: Theory and
  Experiment\/}, \bold{01}, P01005.

Franke, J. et al. (2011). Evolutionary accessibility of mutational
pathways. \emph{PLoS Computational Biology\/}, \bold{7}(8), 1--9.

}

\author{ Ramon Diaz-Uriarte
	NK model implementation by: Ana Hernandez, Alba Álvarez & Coral Fustero.
}

\seealso{
  
  \code{\link{oncoSimulIndiv}},
  \code{\link{plot.genotype_fitness_matrix}},
  \code{\link{evalAllGenotypes}}
  \code{\link{allFitnessEffects}}
  \code{\link{plotFitnessLandscape}}

}
\examples{

## Random fitness for four genes-genotypes,
## plotting and simulating an oncogenetic trajectory

r1 <- rfitness(4)
plot(r1)
oncoSimulIndiv(allFitnessEffects(genotFitness = r1))

#Calculate NK fitness.
nk1 <- rfitness(g = 4, k = 2, k_model = "arithmetic")
plot(nk1)
oncoSimulIndiv(allFitnessEffects(genotFitness = nk1))


}

\keyword{ datagen }

