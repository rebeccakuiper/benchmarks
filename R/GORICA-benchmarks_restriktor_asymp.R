#' Benchmarks for GORIC(A) weights and Posterior Model Probabilities
#'
#' This function calculates case-specific benchmarks for GORIC(A) weights and Posterior Model Probabilities.
#'
#' @param goric_obj 
#' @param pop.est Optional. By default, pop.est = NULL.
#' @param scaling.est Optional. By default, scaling.est = 1. scaling.est = c(.5, 1, 2)
#' @param other.N Optional. By default, other.N = NULL.
#' @param CI Optional. By default, CI = TRUE.
#' @param CI.iter Optional. By default, CI.iter = 1000.
#'
#' @return The output comprises, among other things, the Chi-bar-square weights and critical value for the Chi-bar-square statistic (if asked for) the Chi-bar-square statistic and corresponding p-value.
#' @importFrom restriktor goric
#' @importFrom MASS mvrnorm
# #' @export print.benchmarks
# #' @export summary.benchmarks
#' @export
#' @examples
#'
#' # library(benchmarks)
#'
#' # TO DO Voorbeelden
#' 
#' # TO DO check of dit werkt - ms anova bekijken
#' es_all <- c(0, .2, .5, .8)
#' pop.est <- matrix(NA, nrow = length(es_all), ncol = k)
#' rn <- c("No", "Small", "Medium", "Large")
#' rownames(pop.est) <- rn
#' for(j in 1:length(es_all)){
#'   # j = 1
#'   es_j <- es_all[j]
#'   
#'   i = 1:k
#'   #compute equal difference scores, d, between the means
#'   d =(2*sqrt(k)*es_j) / sqrt(sum((2*i-1-k)^2))
#'   # compute first (read lowest) mean mu
#'   mu.i = (-(k-1)*d) / 2
#'   # compute k means
#'   #means = c(mu.i, mu.i + 1:(k-1)*d) # 1:2:3:...:k-1:k
#'   means_j = c(mu.i + (k-1):1*d, mu.i) # k:k-1:...:3:2:1
#'   
#'   pop.est[j,] <- means_j
#' }
#' colnames(pop.est) <- names(coef(fit))
#'

benchmarks <- function(goric_obj, pop.est = NULL, other.N = NULL, CI.iter = 1000, seed.value = 123) {

  # When testing:
  #goric_obj <- results1
  #pop.es <- 0.2 # c(0, .2) # c(0, .2, .5, .8)
  #ratio.pop.means <- c(3, 2, 1)
  #other.N <- NULL
  # seed.value <- 123
  # CI.iter <- 1000
  
  
  # Check:
  if(!any(class(goric_obj) == "con_goric")){
    return(paste0("The argument goric_obj should be of class con_goric (a goric object from restriktor); it belongs to ", class(goric_obj)))
  }
  
  hypos <- goric_obj$hypotheses_usr
  nr.hypos <- dim(goric_obj$result)[1]
  n.coef <- length(coef(goric_obj))
  nr.es <- length(pop.est)/n.coef
  
  
  if(is.null(pop.est)){
    pop.est <- coef(goric_obj$model.org) # or: goric_obj$model.org$coefficients
  }
  colnames(pop.est) <- names(goric_obj$model.org$coefficients)
  
  #n.coef <- length(b.ratios) # If one of the 2 dim's is 1 - always the case in lm?
  VCOV <- vcov(goric_obj$model.org)
  samplesize <- length(goric_obj$model.org$residuals) # If one of the 2 dim's is 1 - always the case in lm?
  if(!is.null(other.N)){
    VCOV <- VCOV * samplesize / other.N
    samplesize <- other.N
  }
  
  
  nr.iter <- CI.iter
  set.seed(seed.value)
  CI.benchmarks_all <- list(pop.es = pop.es)
  #
  quant <- c(.05, .50, .95)
  CI.benchmarks <- matrix(NA, nrow = nr.hypos, ncol = 1+1+length(quant))
  for(teller.es in 1:nr.es){
    #teller.es = 1
    
    est <- mvrnorm(CI.iter, pop.est[teller.es,], VCOV)
    #
    # Test
    # var(est)
    # VCOV
    # est_test <- mvrnorm(CI.iter*100, pop.est, VCOV)
    # var(est_test)
    # Equals approx VCOV
    
    benchmarks.CI <- matrix(NA, nrow = nr.hypos, ncol = CI.iter)
    for(i in 1:CI.iter){
      # i = 1
      pop.est.CI <- est[i,]
      # Apply GORIC #
      #set.seed(123) 
      results.CI <- goric(pop.est.CI, VCOV = VCOV, constraints = hypos, comparison = "complement", type = "gorica")
      benchmarks.CI[,i] <- matrix(results.CI$result[,7], ncol = 1)
    }
    
    CI.benchmarks_s <- matrix(NA, nrow = nr.hypos, ncol = length(quant))
    Mean <- matrix(NA, nrow = nr.hypos, ncol = 1)
    for(j in 1:nr.hypos){
      # j = 1
      CI.benchmarks_s[j,] <- quantile(benchmarks.CI[j,], quant)
      Mean[j,] <- mean(benchmarks.CI[j,])
    }
    #rownames(CI.benchmarks_s) <- results$result[,1]
    #colnames(CI.benchmarks_s) <- c("5%", "50%", "95%")
    
    CI.benchmarks[,1] <- goric_obj$result[,7] # so in sample
    CI.benchmarks[,2] <- Mean
    CI.benchmarks[,3:(2+length(quant))] <- CI.benchmarks_s
    rownames(CI.benchmarks) <- goric_obj$result[,1]
    colnames(CI.benchmarks) <- c("Sample", "Mean", "5%", "50%", "95%")
    
    name <- paste0("pop.es = ", pop.es[teller.es])
    CI.benchmarks_all[[name]] <- CI.benchmarks
  }
  
  
  benchmarks_all <- matrix(NA, ncol = 1+nr.es, nrow = nr.hypos)
  benchmarks_all[,1] <- CI.benchmarks_all[[2]][,1]
  for(teller.es in 1:nr.es){
    #teller.es = 1
    benchmarks_all[,1+teller.es] <- CI.benchmarks_all[[1+teller.es]][,2]
  }
  colnames(benchmarks_all) <- c("Sample", paste0("pop.es = ", pop.es))
  rownames(benchmarks_all) <- goric_obj$result[,1]
  #benchmarks_all
  
  
  final <- list(#message = message,
    n.coef=k, N = samplesize,
    pop.estimates = pop.est, pop.VCOV = VCOV, 
    benchmarks = benchmarks_all,
    CI.benchmarks = CI.benchmarks_all)
  
  
  class(final) <- c("benchmarks", "list")
  final
  
} # end of function

