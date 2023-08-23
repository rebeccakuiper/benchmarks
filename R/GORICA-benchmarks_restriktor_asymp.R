#' Benchmarks for GORIC(A) weights
#'
#' This function calculates case-specific benchmarks for the GORIC(A) weight and ratio of weights of the preferred hypothesis, assuming user-specified population parameter estimates.
#'
#' @param goric_obj An object from the goric function from the restriktor package. In this function, the GORIC or GORICA can be applied to every type of statistical model.
#' @param pop.est Optional. A nr.es-times-k matrix of population values for the k parameters of interest, leading to nr.es sets of population values. By default, pop.est = NULL; then, the estimates from the sample will be used. Benchmarks will be calculated for each of these value(s).
#' @param other.N Optional (which only works if the goric object is based on a fit object and not on estimates and their covariance matrix). A scalar to denote the (total) sample sizes, if you like to use another than used in the data. You could use this, for instance, to see to which values the GORIC(A) weights will converge (and thus to see the maximum value of the weights). By default, other.N = NULL. In that case, the sample size from the data will be used.
#' @param iter Optional. A scalar denoting the number of iterations used to determine the benchmarks weights. By default, iter = 1000. Notable, the higher iter, the longer the computation time.
#' @param seed.value. Optional. A scalar denoting the seed value. By default, seed.value = 123. By changing this value, you can inspect the sensitivity of the benchmarks -- which should not be present; if it is, you should increase the value of iter.
#'
#' @return The output comprises case-specific benchmarks for the GORIC(A) weight and ratio of weights of the preferred hypothesis (for each of the specified population estimates), based on the 5\%, 35\%, 50\%, 65\%, and 95\% percentiles (using the sample size from the data).
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
#' # Calculate case-specific benchmarks
#' benchmarks_goric_asymp <- benchmarks(goric.obj, pop.est, iter = iter)
#' benchmarks_goric_asymp$error.prob.pref.hypo
#' benchmarks_goric_asymp$pop.estimates
#' benchmarks_goric_asymp$benchmarks.weight
#' benchmarks_goric_asymp$benchmarks.ratios
#'


benchmarks <- function(goric_obj, pop.est = NULL, other.N = NULL, iter = 1000, seed.value = 123) {

  # When testing:
  #goric_obj <- My_goric_obj # results1
  #pop.est <- matrix(c(.22, .51, .83), nrow = 1)
  #other.N <- NULL
  # seed.value <- 123
  # iter <- 1000


  # Check:
  if(!any(class(goric_obj) == "con_goric")){
    return(paste0("The argument goric_obj should be of class con_goric (a goric object from restriktor); it belongs to ", class(goric_obj)))
  }

  hypos <- goric_obj$hypotheses_usr
  nr.hypos <- dim(goric_obj$result)[1]
  PrefHypo <- which.max(goric_obj$result[,7]) #which.max(goric_obj$result$goric.weights)
  pref.hypo <- goric_obj$result$model[PrefHypo]

  n.coef <- length(coef(goric_obj))
  nr.es <- length(pop.est)/n.coef


  #n.coef <- length(b.ratios) # If one of the 2 dim's is 1 - always the case in lm?
  if(is.null(goric_obj$model.org)){
    est_text <- paste0("goric_obj$objectList$", goric_obj$objectNames, "$b.unrestr")
    est_sample <- matrix(eval(parse(text = est_text)), nrow=1)
    if(is.null(pop.est)){
      pop.est <- est_sample
    }
    colnames(pop.est) <- colnames(est_sample)
    #
    #vcov_text <- paste0("goric_obj$objectList$", goric_obj$objectNames, "$CON$VCOV")
    vcov_text <- paste0("goric_obj$objectList$", goric_obj$objectNames, "$Sigma")
    VCOV <- eval(parse(text = vcov_text))
    # Note possible to determine samplesize now (unless part of input), so:
    other.N = NULL
  }else{
    if(is.null(pop.est)){
      pop.est <- coef(goric_obj$model.org) # or: goric_obj$model.org$coefficients
    }
    colnames(pop.est) <- names(goric_obj$model.org$coefficients)
    #
    VCOV <- vcov(goric_obj$model.org)
    samplesize <- length(goric_obj$model.org$residuals) # If one of the 2 dim's is 1 - always the case in lm?
  }
  #
  if(!is.null(other.N)){
    VCOV <- VCOV * samplesize / other.N
    samplesize <- other.N
  }


  nr.iter <- iter
  set.seed(seed.value)
  #
  quant <- c(.05, .35, .50, .65, .95)
  names_quant <- c("Sample", "5%", "35%", "50%", "65%", "95%")
  #
  CI.benchmarks_all <- NULL
  CI.benchmarks_gw_all <- NULL
  for(teller.es in 1:nr.es){
    #teller.es = 1

    est <- mvrnorm(iter, pop.est[teller.es,], VCOV)
    #
    # Test
    # var(est)
    # VCOV
    # est_test <- mvrnorm(iter*100, pop.est, VCOV)
    # var(est_test)
    # Equals approx VCOV

    goric <- rep(NA, nr.iter)
    gw <- matrix(NA, nrow = nr.hypos, ncol = iter)
    for(i in 1:iter){
      # i = 1
      pop.est.CI <- est[i,]
      # Apply GORIC #
      #set.seed(123)
      results.goric <- goric(pop.est.CI, VCOV = VCOV, constraints = hypos, comparison = goric_obj$comparison, type = "gorica")

      goric[i] <- results.goric$result[PrefHypo,7]
      gw[,i] <- results.goric$ratio.gw[PrefHypo,]
    }

    CI.benchmarks_goric <- matrix(c(goric_obj$result[PrefHypo,7], quantile(goric, quant)), nrow = 1) # sample weight with calculated quantiles/percentiles
    colnames(CI.benchmarks_goric) <- names_quant
    rownames(CI.benchmarks_goric) <- pref.hypo
    #
    CI.benchmarks_gw <- matrix(NA, nrow = nr.hypos, ncol = 1+length(quant))
    CI.benchmarks_gw[,1] <- goric_obj$ratio.gw[PrefHypo,] # so in sample
    for(j in 1:nr.hypos){
      CI.benchmarks_gw[j,2:(1+length(quant))] <- quantile(gw[j,], quant)
    }
    colnames(CI.benchmarks_gw) <- names_quant
    rownames(CI.benchmarks_gw) <- paste(pref.hypo, names(goric_obj$ratio.gw[PrefHypo,]))
    #
    #CI.benchmarks_goric
    #CI.benchmarks_gw

    name <- paste0("nr.pop.est = ", teller.es)
    CI.benchmarks_all[[name]] <- CI.benchmarks_goric
    CI.benchmarks_gw_all[[name]] <- CI.benchmarks_gw
    #CI.benchmarks_all
    #CI.benchmarks_gw_all
  }


  # Error probability based on complement of preferred hypothesis in data
  # TO DO zie ANOVA fie
  if(nr.hypos == 2 & goric_obj$comparison == "complement"){
    error.prob <- 1 - goric_obj$result$goric.weights[PrefHypo]
  }else{
    if(PrefHypo == nr.hypos){
      error.prob <- "The failsafe is preferred..."
    }else{
      H_pref <- hypos[[PrefHypo]]
      # Use goric, because ANOVA
      # TO DO what if started with est + cov mx and goric? Dan based on type of input dat gebruiken??
      # Lukt me niet om $b.unrestr en $Sigma te gebruiken...
      #if()
      #
      if(is.null(goric_obj$model.org)){
        results.goric_pref <- goric(est_sample, VCOV = VCOV, constraints = list(H_pref = H_pref), comparison = "complement", type = "goric")
      }else{
        fit_data <- goric_obj$model.org
        results.goric_pref <- goric(fit_data, constraints = list(H_pref = H_pref), comparison = "complement", type = "goric")
      }
      error.prob <- results.goric_pref$result$goric.weights[2]
    }
  }


  final <- list(#message = message,
    n.coef = n.coef, N = samplesize,
    pop.estimates = pop.est, pop.VCOV = VCOV,
    pref.hypo = pref.hypo, error.prob.pref.hypo = error.prob,
    benchmarks.weight = CI.benchmarks_all,
    benchmarks.ratios = CI.benchmarks_gw_all)


  class(final) <- c("benchmarks", "list")
  final

} # end of function

