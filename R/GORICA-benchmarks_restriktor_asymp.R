#' Benchmarks for GORIC(A) weights
#'
#' This function calculates case-specific benchmarks for GORIC(A) weights, together with their 95\% confidence intervals (CIs).
#'
#' @param goric_obj An object from the goric function from the restriktor package. In this function, the GORIC or GORICA can be applied to every type of statistical model.
#' @param pop.est Optional. A nr.es-times-k matrix of population values for the k parameters of interest, leading to nr.es sets of population values. By default, pop.est = NULL; then, the estimates from the sample will be used. Benchmarks will be calculated for each of these value(s).
#' @param other.N Optional. A scalar to denote the (total) sample sizes, if you like to use another than used in the data. You could use this, for instance, to see to which values the GORIC(A) weights will converge (and thus to see the maximum value of the weights). By default, other.N = NULL. In that case, the sample size from the data will be used.
#' @param CI.iter Optional. A scalar denoting the number of iterations used to determine the benchmarks weights and their 95\% confidence intervals (CIs). By default, CI.iter = 1000. Notable, the higher CI.iter, the longer the computation time.
#' @param seed.value. Optional. A scalar denoting the seed value. By default, seed.value = 123. By changing this value, you can inspect the sensitivity of the benchmarks -- which should not be present; if it is, you should increase CI.iter.
#'
#' @return The output comprises case-specific benchmarks for GORIC(A) weights (for each of the specified population estimates) and their 95\% confidence intervals (CIs).
#' @importFrom restriktor goric
#' @importFrom MASS mvrnorm
# #' @export print.benchmarks
# #' @export summary.benchmarks
#' @export
#' @examples
#'
#' # library(benchmarks)
#'
#' # Input 'data'
#' goric.obj <- My_goric_obj
#'
#' # Input w.r.t. population values
#' pop.est <- My_pop_est
#'
#' # Extra
#' iter <- 100
#'
#' # Calculate case-specific benchmarks and their CIs
#' benchmarks_goric_asymp <- benchmarks(My_goric_obj, pop.est, CI.iter = iter)
#' benchmarks_goric_asymp$pop.estimates
#' benchmarks_goric_asymp$benchmarks
#' benchmarks_goric_asymp$CI.benchmarks
#'

benchmarks <- function(goric_obj, pop.est = NULL, other.N = NULL, CI.iter = 1000, seed.value = 123) {

  # When testing:
  #goric_obj <- results1
  #pop.est <-  # c(0, .2, .5, .8)
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

