
prob.k<-function (k, poissonRatePrior) 
{
  Denom <- (poissonRatePrior + 1)^(k + 1)
  Prob <- poissonRatePrior/Denom
  return(Prob)
}

plotPrior_mod<-function (mcmc, expectedNumberOfShifts = 1, burnin = 0.15, priorCol = "red", 
          postCol = "light blue", legendPos = "topright", ...) 
{
  if (!class(mcmc) %in% c("character", "data.frame", "matrix")) {
    stop("mcmc must be either a dataframe or the path to the mcmc_out file.")
  }
  if (is.character(mcmc)) {
    mcmc <- read.csv(mcmc, stringsAsFactors = FALSE)
  }
  mcmc2 <- mcmc[floor(burnin * nrow(mcmc)):nrow(mcmc), ]
  obsK <- seq(from = 0, to = max(mcmc2[, "N_shifts"]), by = 1)
  prior <- sapply(obsK, prob.k, poissonRatePrior = 1/expectedNumberOfShifts)
  prior <- data.frame(N_shifts = obsK, prob = prior)
  posterior <- sapply(obsK, function(x) length(which(mcmc2[, 
                                                           "N_shifts"] == x)))/nrow(mcmc2)
  names(posterior) <- obsK
  posterior <- data.frame(N_shifts = names(posterior), prob = posterior)
  df.posterior<-barplot(posterior[, 2], names.arg = posterior[, 1], ylim = c(0, 
                                                       max(c(prior[, 2], posterior[, 2]))), border = "black", 
          col = postCol,...)
  lines(df.posterior,prior[, 2], col = priorCol)
  points(df.posterior,prior[, 2], col = priorCol, pch=20)
  
  invisible(cbind(N_shifts = prior$N_shifts, priorProbs = prior$prob, 
                  postProbs = posterior$prob))
}

