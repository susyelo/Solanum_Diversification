prob.k<-function (k, poissonRatePrior) 
{
  Denom <- (poissonRatePrior + 1)^(k + 1)
  Prob <- poissonRatePrior/Denom
  return(Prob)
}


computeBayesFactors_mod<-function (postdata, expectedNumberOfShifts, burnin = 0.1, ...) 
{
  if (hasArg("strict") | hasArg("threshpost") | hasArg("threshprior") | 
        hasArg("nbprior") | hasArg("priordata") | hasArg("modelset")) {
    cat("Error - you have specified some argument names that have been deprecated\n")
    cat("in this version of BAMMtools. Check the help file on this function\n")
    cat("to see what has changed\n\n")
    stop()
  }
  if (class(postdata) == "character") {
    dpost <- read.csv(postdata, header = T)
  }
  else if (class(postdata) == "data.frame") {
    dpost <- postdata
  }
  else {
    stop("invalid postdata argument (wrong class) in computeBayesFactors\n")
  }
  dpost <- dpost[floor(burnin * nrow(dpost)):nrow(dpost),]
            
  if (unique(dpost$N_shifts)==0){
    tx <- table(dpost$N_shifts)/nrow(dpost)
    post <- data.frame(N_shifts = as.numeric(names(tx)), prob = as.numeric(tx))
    ux <- as.numeric(names(tx))
    if (length(ux) <= 1) {
      cat("Not enough models sampled in simulation of posterior\n")
      cat("You must have valid posterior probabilities for at least 2 models\n")
      cat("to use this function'\n")
    }
    
  }else
  {
    tx <- table(dpost$N_shifts)/nrow(dpost)
    null_m0<-1/nrow(dpost)
    post <- data.frame(N_shifts = c(0,as.numeric(names(tx))), prob = c(null_m0,as.numeric(tx)))
    ux <- c(0,as.numeric(names(tx)))
    if (length(ux) <= 1) {
      cat("Not enough models sampled in simulation of posterior\n")
      cat("You must have valid posterior probabilities for at least 2 models\n")
      cat("to use this function'\n")
    }
  }
  
  pp <- (1/(1 + expectedNumberOfShifts))
  #prior <- dgeom(ux, prob = pp)
  
  obsK <- seq(from = 0, to = max(dpost[, "N_shifts"]), by = 1)
  prior <- sapply(obsK, prob.k, poissonRatePrior = 1/expectedNumberOfShifts)
  #prior <- data.frame(N_shifts = obsK, prob = prior)
  names(prior) <- obsK
  
  mm <- matrix(NA, nrow = length(prior), ncol = length(prior))
  rownames(mm) <- names(prior)
  colnames(mm) <- names(prior)
  for (i in 1:length(prior)) {
    mi <- ux[i]
    for (j in 1:length(prior)) {
      mj <- ux[j]
      prior_odds <- prior[i]/prior[j]
      post_odds <- post$prob[post$N_shifts == mi]/post$prob[post$N_shifts == 
                                                              mj]
      mm[i, j] <- post_odds * (1/prior_odds)
    }
  }
  return(mm)
}

