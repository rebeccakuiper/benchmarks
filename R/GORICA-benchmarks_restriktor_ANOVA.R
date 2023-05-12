#' Benchmarks for GORIC(A) weights
#'
#' This function calculates case-specific benchmarks for GORIC(A) weights, together with their 95\% confidence intervals (CIs).
#'
#' @param goric_obj An object from the goric function from the restriktor package. In this function, it is assumed that the GORIC is applied to an ANOVA with k group means.
#' @param pop.es Optional. A scalar or vector of population Cohen's d (effect size) value. By default, pop.es = .2. Benchmarks will be calculated for each of these value(s).
#' @param ratio.pop.means Optional. A k-times-1 vector, denoting the relative difference between the k group means. Note that a ratio of c(3,2,1) gives the same as c(1,0,-1), since the consecutive relative differences are 1 in both ratios. every time minus same difference. By default, ratio.pop.means = NULL. In that case, the relative differences from the data are used.
#' @param other.N Optional. A k-times-1 vector or a scalar to denote the group sizes, if you like to use others than used in the data. You could use this, for instance, to see to which values the GORIC(A) weights will converge (and thus to see the maximum value of the weights). If you specify a scalar, it is assumed that each group is of that size. By default, other.N = NULL. In that case, the group sizes from the data will be used.
#' @param CI.iter Optional. A scalar denoting the number of iterations used to determine the benchmarks weights and their 95\% confidence intervals (CIs). By default, CI.iter = 1000. Notable, the higher CI.iter, the longer the computation time.
#' @param seed.value. Optional. A scalar denoting the seed value. By default, seed.value = 123. By changing this value, you can inspect the sensitivity of the benchmarks -- which should not be present; if it is, you should increase CI.iter.
#'
#' @return The output comprises case-specific benchmarks for GORIC(A) weights (for each of the specified effect sizes) and their 95\% confidence intervals (CIs).
#' @importFrom restriktor goric
#' @importFrom fastDummies dummy_cols
# #' @export print.benchmarks_ANOVA
# #' @export summary.benchmarks_ANOVA
#' @export
#' @examples
#'
#' # library(benchmarks)
#'
#' # Input 'data'
#' goric.obj <- My_goric_obj
#'
#' # Input w.r.t. population values
#' pop.es <- c(0, .2, .5, .8)
#' ratio.pop.means <- c(3,2,1)
#'
#' # Extra
#' iter <- 100
#'
#' # Calculate case-specific benchmarks and their CIs
#' benchmarks_goric <- benchmarks_ANOVA(goric.obj, pop.es, ratio.pop.means, CI.iter = iter)
#' benchmarks_goric$pop.means
#' benchmarks_goric$benchmarks
#' benchmarks_goric$CI.benchmarks
#'
#' # An example to see what maximum value of the weights is.
#' # If there is a maximum, then there is overlap in the hypotheses (see guidelines).
#' pop.es <- c(.2, .5, .8)
#' benchmarks_goric_1000 <- benchmarks_ANOVA(goric.obj, pop.es, ratio.pop.means, other.N = 1000, CI.iter = iter)
#' benchmarks_goric_1000$benchmarks
#' benchmarks_goric_1000$CI.benchmarks
#'


benchmarks_ANOVA <- function(goric_obj, pop.es = .2, ratio.pop.means = NULL, other.N = NULL, CI.iter = 1000, seed.value = 123) {

  # When testing:
  #goric_obj <- My_goric_obj # results1
  #pop.es <- 0.2 # c(0, .2) # c(0, .2, .5, .8)
  #ratio.pop.means <- c(3, 2, 1)
  #other.N <- NULL
  # seed.value <- 123
  # CI.iter <- 1000


  # Check:
  if(!any(class(goric_obj) == "con_goric")){
    return(paste0("The argument goric_obj should be of class con_goric (a goric object from restriktor); it belongs to ", class(goric_obj)))
  }
  # TO DO werkt nu niet als class gelijk is aan: "con_gorica" "con_goric"
  # Dan nl model.org niet bekend... lastig met var.e ook


  # number of groups
  n.coef <- length(coef(goric_obj))

  # subject per group
  samplesize <- summary(goric_obj$model.org$model[,2])

  # ES and ratio in data
  # TO DO $b.unrestr en $Sigma zouden het ws moeten zijn, maar die geven dan NULL?
  # Maar in $objectList hebben ze wel waardes....
  means <- coef(goric_obj$model.org) # TO DO model.org bestaat alleen als niet est+vcov input!
  var.e_data <- (sum(goric_obj$model.org$residuals^2) / (sum(samplesize) - n.coef)) # Dit werkt dan ook niet!
  ES_data <- (1/sqrt(var.e_data)) * sqrt((1/n.coef) * sum((means - mean(means))^2))
  #
  ratio_data <- rep(NA, n.coef)
  ratio_data[order(means) == 1] <- 1
  ratio_data[order(means) == 2] <- 2
  d <- means[order(means) == 2] - means[order(means) == 1]
  #d
  #ratio_data
  for(i in 1:n.coef){
    if(order(means)[i] > 2){
      ratio_data[i] <- 1 + (means[i] - means[order(means) == 1])/d
    }
  }
  #ratio_data


  # effect size population
  es <- pop.es
  nr.es <- length(pop.es)

  # ratio population means
  if(is.null(ratio.pop.means)){
    # Then same as in data
    ratio.pop.means <- coef(goric_obj$model.org)
  }else{
    if(length(ratio.pop.means) != n.coef) return(paste0("The argument ratio.pop.means should be of length ", n.coef, " (or NULL) but not of length ", length(ratio.pop.means)))
  }

  # Hypotheses
  hypos <- goric_obj$hypotheses_usr
  nr.hypos <- dim(goric_obj$result)[1]

  # Error variance
  # var.e <- var(goric_obj$model.org$residuals)
  # var.e <- 1
  var.e <- (sum(goric_obj$model.org$residuals^2) / (sum(samplesize) - n.coef))
  #
  # When determining pop. means, value does not matter: works exactly the same
  # choose first or last, then pop. means comparable to sample estimates

  # Possibly adjust var.e based on other sample size
  if(!is.null(other.N)){
    if(length(other.N) == 1){
      var.e <- var.e * (sum(samplesize) - n.coef)
      samplesize <- rep(other.N, n.coef)
      var.e <- var.e / (sum(samplesize) - n.coef)
    }else if(length(other.N) == n.coef){
      var.e <- var.e * (sum(samplesize) - n.coef)
      samplesize <- other.N
      var.e <- var.e / (sum(samplesize) - n.coef)
    }else{
      return(paste0("The argument other.N should be of length 1 or ", n.coef, " (or NULL) but not of length ", length(other.N)))
    }
  }

  means_pop_all <- matrix(NA, ncol = n.coef, nrow = nr.es)
  for(teller.es in 1:nr.es){
    #teller.es = 1

    # Determine mean values, with ratio of ratio.m
    # Solve for x here
    #
    #If all equal, then set population means to all 0
    if(length(unique(ratio.pop.means)) == 1){
      means_pop <- rep(0, n.coef)
    }else{
      fun <- function (d) {
        means_pop = ratio.pop.means*d
        (1/sqrt(var.e)) * sqrt((1/n.coef) * sum((means_pop - mean(means_pop))^2)) - es[teller.es]
      }
      d <- uniroot(fun, lower=0, upper=100)$root
      # Construct means_pop
      means_pop <- ratio.pop.means*d
      # Check
      #means_pop
      #means_pop[1] - means_pop[2]; means_pop[2] - means_pop[3]
      # Check: ES =
      #(1/sqrt(var.e)) * sqrt((1/n.coef) * sum((means_pop - mean(means_pop))^2))
    }

    means_pop_all[teller.es, ] <- means_pop
  }
  colnames(means_pop_all) <- colnames(coef(goric_obj))
  rownames(means_pop_all) <- paste0("pop.es = ", pop.es)
  # means_pop_all

  # Create dummies
  sample <- NULL
  D_ <- matrix(NA, nrow = sum(samplesize))
  D_[1:samplesize[1]] <- 1
  for(j in 2:length(samplesize)){
    # j = 2
    D_[(1+sum(samplesize[1:(j-1)])):sum(samplesize[1:j])] <- j
  }
  # D_
  D <- as.factor(D_)
  sample$D <- D
  sample <- dummy_cols(sample, select_columns = 'D')



  nr.iter <- CI.iter

  set.seed(seed.value)
  CI.benchmarks_all <- list(pop.es = pop.es)
  for(teller.es in 1:nr.es){
    #teller.es = 1
    means_pop <- means_pop_all[teller.es, ]

    #means_est <- matrix(NA, ncol = n.coef, nrow = nr.iter)
    #covmx_est <- array(NA, dim = c(n.coef, n.coef, nr.iter))
    goric <- matrix(NA, ncol = nr.hypos, nrow = nr.iter)
    #gorica <- matrix(NA, ncol = nr.hypos, nrow = nr.iter)
    for(i in 1:nr.iter){
      # teller.es = 1; i = 1

      # Sample residuals
      epsilon <- rnorm(sum(samplesize), sd=sqrt(var.e))
      #
      # Generate data
      sample$y <- as.matrix(sample[,2:(1+n.coef)]) %*% matrix(means_pop, nrow = n.coef) + epsilon

      # Obtain fit
      fit <- lm(y ~ 0 + D, data=sample)
      #fit


      # Store output

      #means_est[i,] <- coef(fit)
      #covmx_est[, ,i] <- vcov(fit) # Only influenced by res var in data, not means_pop
      #
      # check res. var. in data
      # (sum(fit$residuals^2) / (sum(samplesize) - n.coef)) / samplesize
      # vcov(fit)
      # The same values indeed
      #
      # Btw1: var(epsilon) differs a bit from (sum(fit$residuals^2) / (sum(samplesize) - n.coef))
      #
      ## Btw2: Not: sum(fit$residuals^2) / (sum(samplesize) - n.coef - 1)

      #set.seed(123)
      # GORICA or GORICA depending on what is done in data
      results.goric <- goric(fit, constraints = hypos, comparison = goric_obj$comparison, type = goric_obj$type)
      # GORIC
      #results.goric <- goric(fit, constraints = hypos, comparison = goric_obj$comparison, type = "goric")
      # GORICA
      #results.gorica <- goric(fit, constraints = hypos, comparison = goric_obj$comparison, type = "gorica")
      #
      #Test
      #goric(fit, constraints = hypos, comparison = goric_obj$comparison, type = "gorica")
      #goric(coef(fit), VCOV = vcov(fit), constraints = hypos, comparison = goric_obj$comparison, type = "gorica")

      goric[i,] <- results.goric$result[,7]
      #gorica[i,] <- results.goric$result[,7]
    }

    quant <- c(.05, .50, .95)
    CI.benchmarks_goric <- matrix(NA, nrow = nr.hypos, ncol = 1+1+length(quant))
    #CI.benchmarks_gorica <- matrix(NA, nrow = nr.hypos, ncol = 1+1+length(quant))
    CI.benchmarks_goric[,1] <- goric_obj$result[,7] # so in sample
    #CI.benchmarks_gorica[,1] <- goric_obj$result[,7] # NB dan geen gorica
    CI.benchmarks_goric[,2] <- colMeans(goric)
    #CI.benchmarks_gorica[,2] <- colMeans(gorica)
    for(j in 1:nr.hypos){
      CI.benchmarks_goric[j,3:(2+length(quant))] <- quantile(goric[,j], quant)
      #CI.benchmarks_gorica[j,3:(2+length(quant))] <- quantile(goric[,j], quant)
    }
    colnames(CI.benchmarks_goric) <- c("Sample", "Mean", "5%", "50%", "95%")
    #colnames(CI.benchmarks_gorica) <- c("Sample", "Mean", "5%", "50%", "95%")
    rownames(CI.benchmarks_goric) <- results.goric$result[,1]
    #rownames(CI.benchmarks_gorica) <- results.gorica$result[,1]
    #CI.benchmarks_goric
    #CI.benchmarks_gorica

    name <- paste0("pop.es = ", pop.es[teller.es])
    CI.benchmarks_all[[name]] <- CI.benchmarks_goric
    #CI.benchmarks_all

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


  # Error probability based on complement of preferred hypothesis in data
  PrefHypo <- which.max(goric_obj$result[,7]) #which.max(goric_obj$result$goric.weights)
  pref.hypo <- goric_obj$result$model[PrefHypo]
  if(PrefHypo > nr.hypos){
    error.prob <- "The failsafe is preferred..."
  }else{
    H_pref <- hypos[[PrefHypo]]
    # Use goric, because ANOVA
    # TO DO what if started with est + cov mx and goric? Dan based on type of input dat gebruiken??
    # Lukt me niet om $b.unrestr en $Sigma te gebruiken...
    #if()
    fit_data <- goric_obj$model.org
    results.goric_pref <- goric(fit_data, constraints = list(H_pref = H_pref), comparison = "complement", type = "goric")
    error.prob <- results.goric_pref$result$goric.weights[2]
  }



  final <- list(#message = message,
                n.coef = n.coef,
                group.size = samplesize,
                means.data = means, ratio.means.data = ratio_data, ES.data = ES_data,
                res.var.data = var.e_data,
                pop.es = pop.es, pop.means = means_pop_all,
                ratio.pop.means = ratio.pop.means,
                res.var.pop = var.e,
                pref.hypo = pref.hypo, error.prob.pref.hypo = error.prob,
                benchmarks = benchmarks_all,
                CI.benchmarks = CI.benchmarks_all)

  class(final) <- c("benchmarks", "list")
  final

} # end of function

