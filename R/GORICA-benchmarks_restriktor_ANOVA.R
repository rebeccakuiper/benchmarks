#' Benchmarks for GORIC(A) weights
#'
#' This function calculates, for an ANOVA model, case-specific benchmarks for the GORIC(A) weight and ratio of weights of the preferred hypothesis, assuming user-specified population parameter estimates.
#'
#' @param goric_obj An object from the goric function from the restriktor package. In this function, it is assumed that the GORIC is applied to an ANOVA with k group means.
#' @param pop.es Optional. A scalar or vector of population Cohen's f (effect size) value. By default, pop.es = 0. Benchmarks will be calculated for each of these value(s).
#' @param ratio.pop.means Optional. A k-times-1 vector, denoting the relative difference between the k group means. Note that a ratio of c(3,2,1) gives the same as c(1,0,-1), since the consecutive relative differences are 1 in both ratios. every time minus same difference. By default, ratio.pop.means = NULL. In that case, the relative differences from the data are used.
#' @param N Needed (and only used) if goric object is based on estimates and their covariance matrix (instead of on a model / fit object). A k-times-1 vector or a scalar to denote the group sizes. If you specify a scalar, it is assumed that each group is of that size. By default, N = NULL.
#' @param other.N Optional. A k-times-1 vector or a scalar to denote the group sizes, if you like to use others than used in the data. You could use this, for instance, to see to which values the GORIC(A) weights will converge (and thus to see the maximum value of the weights). If you specify a scalar, it is assumed that each group is of that size. By default, other.N = NULL. In that case, the group sizes from the data will be used.
#' @param iter Optional. A scalar denoting the number of iterations used to determine the benchmarks weights. By default, iter = 1000. Notable, the higher iter, the longer the computation time.
#' @param seed.value. Optional. A scalar denoting the seed value. By default, seed.value = 123. By changing this value, you can inspect the sensitivity of the benchmarks -- which should not be present; if it is, you should increase the value of iter.
#'
#' @return The output comprises case-specific benchmarks for the GORIC(A) weight and ratio of weights of the preferred hypothesis (for each of the specified population estimates), based on the 5\%, 35\%, 50\%, 65\%, and 95\% percentiles (using the sample size from the data).
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
#' pop.es <- c(0, .1, .25, .4) # According to Cohen 1992
#' ratio.pop.means <- c(3,2,1)
#'
#' # Extra
#' iter <- 100
#'
#' # Calculate case-specific benchmarks
#' benchmarks_goric <- benchmarks_ANOVA(goric.obj, pop.es, ratio.pop.means, iter = iter)
#' benchmarks_goric$pop.means
#' benchmarks_goric$benchmarks.weight
#' benchmarks_goric$benchmarks.ratios
#'
#' # An example to see what maximum value of the weights is.
#' # If there is a maximum, then there is overlap in the hypotheses (see guidelines).
#' pop.es <- c(.1, .25, .4)
#' benchmarks_goric_1000 <- benchmarks_ANOVA(goric.obj, pop.es, ratio.pop.means, other.N = 1000, iter = iter)
#' benchmarks_goric_1000$benchmarks.weight
#' benchmarks_goric_1000$benchmarks.ratios
#'


benchmarks_ANOVA <- function(goric_obj, pop.es = 0, ratio.pop.means = NULL, N = NULL, other.N = NULL, iter = 1000, seed.value = 123) {

  # When testing:
  #
  #goric_obj <- My_goric_obj # goric_obj <- results1
  # goric_obj <- output_gorica_c; N <- 30 # N <- rep(30,5)
  # goric_obj <- output_gorica_c_fit; N <- NULL
  # goric_obj <- results_1c
  #
  #pop.es <- 0.2 # pop.es <- c(0, .1) # pop.es <- c(0, .1, .25, .4)
  #ratio.pop.means <- c(3, 2, 1) # ratio.pop.means <- NULL
  #N <- NULL
  #other.N <- NULL
  # seed.value <- 123
  # iter <- 3 # iter <- 100
  # library(fastDummies)


  # Check:
  if(!any(class(goric_obj) == "con_goric")){
    return(paste0("The argument goric_obj should be of class con_goric (a goric object from restriktor); it belongs to ", class(goric_obj)))
  }


  # number of groups
  n.coef <- length(coef(goric_obj))


  # ES and ratio in data
  #
  if(is.null(goric_obj$model.org)){
    # Number of subjects per group
    if(is.null(N)){
      other.N = NULL
      samplesize <- "Could not be retrieved from the input."
      # TO DO MAAK ERROR
    }else if(length(N) == 1){
      samplesize <- rep(N, n.coef)
    }else{
      samplesize <- N
    }

    # Unrestricted means
    est_text <- paste0("goric_obj$objectList$", goric_obj$objectNames, "$b.unrestr")
    means <- eval(parse(text = est_text))

    # residual variance
    vcov_text <- paste0("goric_obj$objectList$", goric_obj$objectNames, "$Sigma")
    VCOV <- eval(parse(text = vcov_text))
    var.e_data_mx <- VCOV * samplesize
    var.e_data <- var.e_data_mx[1,1] # TO DO check always same elements on diagonal?

  }else{
    # Number of subjects per group
    samplesize <- summary(goric_obj$model.org$model[,2])

    # Unrestricted means
    means <- coef(goric_obj$model.org)

    # residual variance
    var.e_data <- (sum(goric_obj$model.org$residuals^2) / (sum(samplesize) - n.coef))
  }


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
    ratio.pop.means <- ratio_data #coef(goric_obj$model.org)
  }else{
    if(length(ratio.pop.means) != n.coef) return(paste0("The argument ratio.pop.means should be of length ", n.coef, " (or NULL) but not of length ", length(ratio.pop.means)))
  }

  # Hypotheses
  hypos <- goric_obj$hypotheses_usr
  nr.hypos <- dim(goric_obj$result)[1]
  PrefHypo <- which.max(goric_obj$result[,7]) #which.max(goric_obj$result$goric.weights)
  pref.hypo <- goric_obj$result$model[PrefHypo]

  # Error variance
  # var.e <- var(goric_obj$model.org$residuals)
  # var.e <- 1
  var.e <- var.e_data
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
  # D_ <- matrix(NA, nrow = sum(samplesize))
  # D_[1:samplesize[1]] <- 1
  # for(j in 2:length(samplesize)){
  #   # j = 2
  #   D_[(1+sum(samplesize[1:(j-1)])):sum(samplesize[1:j])] <- j
  # }
  # # D_
  # D <- as.factor(D_)
  D <- as.factor(rep(1:n.coef, times = samplesize))
  sample$D <- D
  sample <- dummy_cols(sample, select_columns = 'D')
  colnames(sample)[-1] <- names(coef(goric_obj))


  nr.iter <- iter
  set.seed(seed.value)
  #
  quant <- c(.05, .35, .50, .65, .95)
  names_quant <- c("Sample", "5%", "35%", "50%", "65%", "95%")
  #
  CI.benchmarks_all <- NULL
  CI.benchmarks_gw_all <- NULL
  CI.benchmarks_lw_all <- NULL
  CI.benchmarks_lw_ge1_all <- NULL
  CI.benchmarks_ld_all <- NULL
  CI.benchmarks_ld_ge0_all <- NULL
  for(teller.es in 1:nr.es){
    #teller.es = 1
    means_pop <- means_pop_all[teller.es, ]

    #means_est <- matrix(NA, ncol = n.coef, nrow = nr.iter)
    #covmx_est <- array(NA, dim = c(n.coef, n.coef, nr.iter))
    goric <- rep(NA, nr.iter)
    gw <- matrix(NA, nrow = nr.hypos, ncol = iter)
    lw <- matrix(NA, nrow = nr.hypos, ncol = iter)
    ld <- matrix(NA, nrow = nr.hypos, ncol = iter)
    for(i in 1:nr.iter){
      # teller.es = 1; i = 1

      # Sample residuals
      epsilon <- rnorm(sum(samplesize), sd=sqrt(var.e))
      #
      # Generate data
      sample$y <- as.matrix(sample[,2:(1+n.coef)]) %*% matrix(means_pop, nrow = n.coef) + epsilon
      df <- data.frame(y = sample$y, sample[,2:(1+n.coef)])

      # Obtain fit
      #fit <- NULL
      fit <- lm(y ~ 0 + ., data = df)
      #fit
      #summary(fit)


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
      results.goric <- goric(fit,
                             hypotheses = hypos,
                             comparison = goric_obj$comparison,
                             type = goric_obj$type)
      #
      #Test
      #goric(fit, hypotheses = hypos, comparison = goric_obj$comparison, type = "gorica")
      #goric(coef(fit), VCOV = vcov(fit), hypotheses = hypos, comparison = goric_obj$comparison, type = "gorica")

      goric[i] <- results.goric$result[PrefHypo,7]
      gw[,i] <- results.goric$ratio.gw[PrefHypo,]
      lw[,i] <- results.goric$ratio.lw[PrefHypo,]
      ld[,i] <- (results.goric$result$loglik[PrefHypo] - results.goric$result$loglik)
    }

    CI.benchmarks_goric <- matrix(c(goric_obj$result[PrefHypo,7], quantile(goric, quant)), nrow = 1) # sample weight with calculated quantiles/percentiles
    colnames(CI.benchmarks_goric) <- names_quant
    rownames(CI.benchmarks_goric) <- pref.hypo
    #
    #
    CI.benchmarks_gw <- matrix(NA, nrow = nr.hypos, ncol = 1+length(quant))
    CI.benchmarks_lw <- matrix(NA, nrow = nr.hypos, ncol = 1+length(quant))
    CI.benchmarks_lw_ge1 <- matrix(NA, nrow = nr.hypos, ncol = 1+length(quant))
    CI.benchmarks_ld <- matrix(NA, nrow = nr.hypos, ncol = 1+length(quant))
    CI.benchmarks_ld_ge0 <- matrix(NA, nrow = nr.hypos, ncol = 1+length(quant))
    #
    CI.benchmarks_gw[,1] <- goric_obj$ratio.gw[PrefHypo,] # so in sample
    CI.benchmarks_lw[,1] <- goric_obj$ratio.lw[PrefHypo,] # so in sample
    for(j in 1:nr.hypos){
      if(goric_obj$ratio.lw[PrefHypo,j] >= 1){
        CI.benchmarks_lw_ge1[j,1] <- goric_obj$ratio.lw[PrefHypo,j] # so in sample
      }else{
        CI.benchmarks_lw_ge1[j,1] <- 1/goric_obj$ratio.lw[PrefHypo,j] # so in sample
      }
    }
    CI.benchmarks_ld[,1] <- (goric_obj$result$loglik[PrefHypo] - goric_obj$result$loglik) # so in sample
    CI.benchmarks_ld_ge0[,1] <- abs(goric_obj$result$loglik[PrefHypo] - goric_obj$result$loglik) # so in sample
    #
    lw_ge1 <- lw
    lw_ge1[lw < 1] <- 1/lw[lw < 1]
    ld_ge0 <- abs(ld)
    for(j in 1:nr.hypos){
      CI.benchmarks_gw[j,2:(1+length(quant))] <- quantile(gw[j,], quant)
      CI.benchmarks_lw[j,2:(1+length(quant))] <- quantile(lw[j,], quant)
      CI.benchmarks_lw_ge1[j,2:(1+length(quant))] <- quantile(lw_ge1[j,], quant)
      CI.benchmarks_ld[j,2:(1+length(quant))] <- quantile(ld[j,], quant)
      CI.benchmarks_ld_ge0[j,2:(1+length(quant))] <- quantile(ld_ge0[j,], quant)
    }
    #
    colnames(CI.benchmarks_gw) <- colnames(CI.benchmarks_lw) <- colnames(CI.benchmarks_lw_ge1) <- colnames(CI.benchmarks_ld) <- colnames(CI.benchmarks_ld_ge0) <- names_quant
    #
    rownames(CI.benchmarks_gw) <- rownames(CI.benchmarks_lw) <- rownames(CI.benchmarks_lw_ge1) <- rownames(CI.benchmarks_ld) <- rownames(CI.benchmarks_ld_ge0) <- paste(pref.hypo, names(goric_obj$ratio.gw[PrefHypo,]))
    #
    #CI.benchmarks_goric
    #CI.benchmarks_gw
    #CI.benchmarks_lw
    #CI.benchmarks_lw_ge1
    #CI.benchmarks_ld
    #CI.benchmarks_ld_ge0

    name <- paste0("pop.es = ", pop.es[teller.es])
    CI.benchmarks_all[[name]] <- CI.benchmarks_goric
    CI.benchmarks_gw_all[[name]] <- CI.benchmarks_gw
    CI.benchmarks_lw_all[[name]] <- CI.benchmarks_lw
    CI.benchmarks_lw_ge1_all[[name]] <- CI.benchmarks_lw_ge1
    CI.benchmarks_ld_all[[name]] <- CI.benchmarks_ld
    CI.benchmarks_ld_ge0_all[[name]] <- CI.benchmarks_ld_ge0
    #CI.benchmarks_all
    #CI.benchmarks_gw_all
    #CI.benchmarks_lw_all
    #CI.benchmarks_lw_ge1_all
    #CI.benchmarks_ld_all
    #CI.benchmarks_ld_ge0_all

  }


  # Error probability based on complement of preferred hypothesis in data

  if(nr.hypos == 2 & goric_obj$comparison == "complement"){
    if(goric_obj$type == 'goric'){
      error.prob <- 1 - goric_obj$result$goric.weights[PrefHypo]
    }else{
      error.prob <- 1 - goric_obj$result$gorica.weights[PrefHypo]
    }
  }else{
    if(PrefHypo == nr.hypos & goric_obj$comparison == "unconstrained"){
      error.prob <- "The unconstrained (i.e., the failsafe) containing all possible orderings is preferred..."
    }else{
      H_pref <- hypos[[PrefHypo]]
      if(is.null(goric_obj$model.org)){
        results.goric_pref <- goric(means, VCOV = VCOV,
                                    hypotheses = list(H_pref = H_pref),
                                    comparison = "complement",
                                    type = goric_obj$type)
      }else{
        fit_data <- goric_obj$model.org
        results.goric_pref <- goric(fit_data,
                                    hypotheses = list(H_pref = H_pref),
                                    comparison = "complement",
                                    type = goric_obj$type)
        # Use same type as in original model (could do goric)
      }
      if(goric_obj$type == 'goric'){
        error.prob <- results.goric_pref$result$goric.weights[2]
      }else{
        error.prob <- results.goric_pref$result$gorica.weights[2]
      }
    }
  }
  #
  # TO DO bepaal ook quantiles voor error prob. Ws verwerken in bovenstaande!
  # Is het zinnig? Op zich zou error prob al zinnig moeten zijn immers....



  final <- list(#message = message,
                n.coef = n.coef,
                group.size = samplesize,
                means.data = means, ratio.means.data = ratio_data, ES.data = ES_data,
                res.var.data = var.e_data,
                pop.es = pop.es, pop.means = means_pop_all,
                ratio.pop.means = ratio.pop.means,
                res.var.pop = var.e,
                pref.hypo = pref.hypo, error.prob.pref.hypo = error.prob,
                benchmarks.weight = CI.benchmarks_all,
                benchmarks.ratios = CI.benchmarks_gw_all,
                benchmarks.LLratios = CI.benchmarks_lw_all,
                benchmarks.LLratios_ge1 = CI.benchmarks_lw_ge1_all,
                benchmarks.difLL = CI.benchmarks_ld_all,
                benchmarks.absdifLL = CI.benchmarks_ld_ge0_all)

  class(final) <- c("benchmarks", "list")
  final

} # end of function

