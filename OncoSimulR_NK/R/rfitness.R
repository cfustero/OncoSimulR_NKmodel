## Copyright 2013, 2014, 2015, 2016, 2017 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


rfitness <- function(g, c= 0.5,
                     sd = 1,
                     mu = 1,
                     reference = "random", ## "random", "max", or the vector,
                                     ## e.g., rep(1, g). If random, a
                                     ## random genotype is chosen as
                                     ## reference. If "max" this is rep(1, g)
                     k = NULL, ## neihgbours
                     k_model = c(NA, "arithmetic", "geometric"), ## Nk algorithm options
                     scale = NULL, ## a two-element vector: min and max
                     wt_is_1 = c("subtract", "divide", "force", "no"),
                     ## wt_is_1 = TRUE, ## wt has fitness 1
                     log = FALSE, ## log only makes sense if all values >
                                 ## 0. scale with min > 0, and/or set
                                 ## wt_is_1 = divide
                     min_accessible_genotypes = NULL,
                     accessible_th = 0,
                     truncate_at_0 = TRUE) {
    ## Like Franke et al., 2011 and others of Krug. Very similar to Greene
    ## and Crona, 2014. And this allows moving from HoC to purely additive
    ## changing c and sd.

    ## FIXME future: do this with order too?
    ##    - do not generate a matrix of genotypes but a matrix of ordered genot.
    wt_is_1 = match.arg(wt_is_1)
    
    m <- generate_matrix_genotypes(g)
    done <- FALSE
    ## attempts <- 0 ## for debugging/tracking purposes
    while(!done) {
        ## attempts <- attempts + 1
      if(is.null(k)){
        f_r <- rnorm(nrow(m), mean = mu, sd = sd)
        
        if(inherits(reference, "character") && length(reference) == 1) {
            if(reference == "random") {
                referenceI <- m[sample(nrow(m), 1), ]
            } else if(reference == "max") {
                referenceI <- rep(1, g)
            } else if(reference == "random2") {
                referenceI <- create_eq_ref(g)
            }
        } else {
            referenceI <- reference
            }
        d_reference <- apply(m, 1, function(x) sum(abs(x - referenceI)))
        f_det <- -c * d_reference
        ## f_det <- rowSums(m) * slope/nrow(m) ## this is Greene and Krona
        fi <- f_r + f_det
      } else{
          ## NK algorithm

    #Firstly we generate the genotypes matrix.
 	  #Now, we create a new matrix, only includying the possible interactions between the locus and its neighbours.
 	  #This means, if K = 2, then the matrix will have three positions and all its possible 0/1 combinations.
       vector_test = NULL
       for (i in 1:length(colnames(m))) {
         vector_test <- c(vector_test, sample(colnames(m)[-i], k))
       }
      
 	  #arrange neighbour matrix
       m_N = matrix(vector_test, nrow=g, byrow=TRUE)
       m_N = as.data.frame(m_N)
       m_N[,(k+1)] = colnames(m)
       m_N= t(apply(m_N,1,sort))
      
       #create reference matrix of genotypes without fitness value
       m_map =  as.data.frame(generate_matrix_genotypes(k+1))
       m_map$fitness <- rep(0,nrow(m_map))
      
       #declare fitness matrix
       #stores finall fitness values
       m_rfit <- matrix(nrow=nrow(m), ncol=ncol(m))
      
      
       #for loop iterate through each row in neighbour matrix: each genotype
         #associate genotype-fitness to depending on the chosen k_model reference_matrix  
       #first apply iterates through rows of fitness matrix 
       #second apply iterates through rows of reference matrix, if reference genotypes matchs with 
         #real situation of genotypes, value is stored
       for (it in 1:dim(m_N)[1]){
         ifelse(k_model == "arithmetic",
                m_map$fitness <- rnorm(nrow(m_map),mean = 0, sd = 1),
                m_map$fitness <- runif(nrow(m_map),min = 0.05, max = 1))
         colnames(m_map) <- c(m_N[it,],"fitness")
         m_rfit[,it] = unlist(apply(m[,c(m_N[it,])],1, 
                                    function(x) apply(m_map, 1, 
                                                      function(z) if(all(x==z[c(1:(k+1))])==TRUE)
                                                      {z[(k+2)]}) ))
       }
      
       #calculates fitness depending on chosen model(geometric or additive)
       if(k_model == "arithmetic") {
         #m = cbind(m,  gFit= apply(m_rfit,1, sum)) 
         fi <- apply(m_rfit,1, sum)
       } else if (k_model == "geometric") {
         #m = cbind(m,  gFit= apply(m_rfit,1, function(x) prod(x)^(1/g)))
         fi <- apply(m_rfit,1, function(x) prod(x)^(1/g))
       }  else {
         warning(' model ', k_model,' is not suitable for NK algorithm' )
       }
      
       #fi <- m[,(g+1)]	
      }
       
        if(!is.null(scale)) {
            fi <- (fi - min(fi))/(max(fi) - min(fi))
            fi <- scale[1] + fi * (scale[2] - scale[1])
        }
      
        if(wt_is_1 == "divide") {
            ## We need to shift to avoid ratios involving negative numbers and
            ## we need to avoid having any fitness at 0, specially the wt.  If
            ## any negative value, add the min, and shift by the magnitude of
            ## the min to avoid any at 0.

            ## If you use scale and wt_is_1, this will move the scale. It is
            ## not possible to obtain a linear transformation that will keep
            ## the min and max of the scale, and wt at 1.
            min_fi <- min(fi)
            if(min_fi < 0)
                fi <- fi + 2 * abs(min(fi))
            fi <- fi/fi[1]
        } else if (wt_is_1 == "subtract") {
            fi <- fi - fi[1] + 1.0
        } else if(wt_is_1 == "force") {
            fi[1] <- 1.0
            if(!is.null(scale)) {
                if( (1 > scale[2]) || (1 < scale[1]))
                    warning("Using wt_is_1 = force and scale, but scale does ",
                            "not include 1")
            }
        }
      
        if(truncate_at_0) {
            fi[fi <= 0] <- 1e-9
        }
        if(log) {
            fi <- log(fi/fi[1]) + 1
        }
        m <- cbind(m, Fitness = fi)
        if(!is_null_na(min_accessible_genotypes)) {
            ## num_accessible_genotypes <- count_accessible_g(m, accessible_th)
            ## Recall accessibleGenotypes includes the wt, if accessible.
            num_accessible_genotypes <- length(wrap_accessibleGenotypes(m, accessible_th)) - 1
            ## message("\n     num_accessible_genotypes = ", num_accessible_genotypes, "\n")
            if(num_accessible_genotypes >= min_accessible_genotypes) {
                done <- TRUE
                attributes(m) <- c(attributes(m),
                                   accessible_genotypes = num_accessible_genotypes,
                                   accessible_th = accessible_th)
            } else {
                ## Cannot start again with a fitness column
                m <- m[, -ncol(m), drop = FALSE]
            }
        } else {
            done <- TRUE
        }
    }
    ## message("\n number of attempts = ", attempts, "\n")
    class(m) <- c(class(m), "genotype_fitness_matrix")
    return(m)
}


create_eq_ref <- function(g) {
    ## "random" gives more prob. to genotypes with
    ## number of mutated genes close to g/2.
    ## This gives equal prob to having the reference
    ## be of any of the possible number of mutated genes.
    nm <- sample(g, 1)
    ref <- c(rep(1, nm), rep(0, g - nm))
    sample(ref)
}

