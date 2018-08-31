
#' @examples
#' ex <- example_finemap()
#' se <- rep(1/sqrt(ex$n1), nrow(ex$tab1))
#' effect <- ex$tab1$zscore * se
#' out <- abf(effect, se)
#' @export
abf <- function(effect, se)
{
  # copmute ABF
  W <- 0.21^2 # suggested by Wakefield
  pi1 <- 1 / length(effect) # prior on alternative (association)

  tab <- BFDPfunV(effect, se^2, W, pi1) %>% bind_cols

  ### add input to `tab`
  ss <- data_frame(effect = effect, se = se) %>%
    mutate(zscore = effect / se)
  
  tab <- bind_cols(ss, tab)

  tab <- mutate(tab, 
    ABF = ABFfun(zscore, se^2, W),
    ABFinv = (1 / ABF))

  tab <- within(tab, sumABFinv <- sum(ABFinv))

  tab <- mutate(tab, P1_ABF = ABFinv / sumABFinv)

  # sort by `P1_ABF`
  tab <- arrange(tab, desc(P1_ABF))
  
  return(tab)
}


#------------------
# ABF computation
#------------------

# This function calculates BFDP, the approximate Pr( H0 | thetahat )
#
# http://faculty.washington.edu/jonno/BFDP.R
#
# doi.org/10.1086/519024, p. 120, formula (6) and above
#
# @param thetahat an estiamte of the log relative risk
# @param V the variance of this estimate
# @param W the prior variance
# @param pi1 the prior probability of a non-null association
BFDPfunV <- function(thetahat, V, W, pi1)
{
  pH0 <- dnorm(thetahat, m = 0, s = sqrt(V))
  postvar <- V + W
  pH1 <- dnorm(thetahat, m = 0, s = sqrt(postvar))
  BF <- pH0/pH1
  
  PO <- (1-pi1) / pi1
  BFDP <- BF * PO / (BF*PO + 1)
  
  list(BF = BF, pH0 = pH0, pH1 = pH1, BFDP = BFDP)
}

# This function computes ABF according to Wakefield (doi.org/10.1086/519024)
#
# doi.org/10.1038/ng.2435, Online Methods, Section Bayesian approach
ABFfun <- function(Z, V, W)
{
  r <- W / (V + W)
  
  (1 / sqrt(1 - r)) * exp(-Z*Z*r/2)
}
