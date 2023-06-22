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
#' @importFrom fastDummies dummy_cols
# #' @export print.benchmarks
# #' @export summary.benchmarks
#' @export
#' @examples
#'
#' # library(benchmarks)
#'
#' # TO DO Voorbeelden
#'

# TO DO 
# ratio c(3,2,1) gives same as c(1,0,-1) -- every time minus same difference
# other.N als 1 getal, dan dat getal per groep; anders juist length n.coef

benchmarks_ANOVA <- function(goric_obj, pop.es = .2, ratio.pop.means = NULL, other.N = NULL, CI.iter = 1000, seed.value = 123) {

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

  
  # effect size
  es <- pop.es 
  nr.es <- length(pop.es)
  # number of groups
  n.coef <- length(coef(goric_obj))
  
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

  # subject per group
  samplesize <- summary(goric_obj$model.org$model[,2])
  
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
      results.goric <- goric(fit, constraints = hypos, comparison = "complement", type = "goric")
      #results.gorica <- goric(fit, constraints = hypos, comparison = "complement", type = "gorica")
      #
      #Test
      #goric(fit, constraints = hypos, comparison = "complement", type = "gorica")
      #goric(coef(fit), VCOV = vcov(fit), constraints = hypos, comparison = "complement", type = "gorica")
          
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

  
  
  final <- list(#message = message,
                n.coef = n.coef,
                group.size = samplesize,
                pop.es = pop.es, pop.means = means_pop_all,
                ratio.pop.means = ratio.pop.means, res.var = var.e,
                benchmarks = benchmarks_all,
                CI.benchmarks = CI.benchmarks_all)
  
  class(final) <- c("benchmarks", "list")
  final
  
} # end of function

