---
title: "Fine-mapping using Approximate Bayes Factor (ABF) on toy example"
author: "Andrey Ziyatdinov"
date: "2018-08-17"
output:
  html_document:
    theme: united
    toc: true
    number_sections: false
    keep_md: true
---



# About

Load libraries:


```r
library(magrittr)
library(dplyr)  

library(ggplot2)
theme_set(theme_minimal())

library(devtools)
load_all("~/git/variani/finemapr/")
```

# Wakefield's functions


```r
# This function calculates BFDP, the approximate Pr( H0 | thetahat )
#
# @href http://faculty.washington.edu/jonno/BFDP.R
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

ABFfun <- function(Z, V, W)
{
  r <- W / (V + W)
  
  (1 / sqrt(1 - r)) * exp(-Z*Z*r/2)
}
```

# Example data set


```r
ex <- example_finemap()

ss <- ex$tab1
ld <- ex$ld1 # not needed for ABF, be used later on for FINEMAP
N <- 4444

# simulate standard errors (se) and add effect sizes (effect)
set.seed(1)
ss <- mutate(ss,
  se = rnorm(n(), 1/sqrt(N), 0.005),
  effect = zscore * se)

# remove strong signals
ss <- filter(ss, abs(zscore) < 6)

# sort
ss <- arrange(ss, -abs(zscore)) %>%
  mutate(rank = seq(1, n())) %>% select(rank, everything())
```

# ABF


```r
W <- 0.21^2 # suggested by Wakefield
pi1 <- 1 / nrow(ss) # prior on alternative (association)

abf <- BFDPfunV(ss$effect, ss$se^2, W, pi1) %>% bind_cols

abf <- bind_cols(ss, abf)

abf <- mutate(abf, 
  ABF = ABFfun(zscore, se^2, W),
  ABFinv = (1 / ABF))
abf <- within(abf, sumABFinv <- sum(ABFinv))

abf <- mutate(abf, P1_ABF = ABFinv / sumABFinv)
```

Top 5 SNPs:


```r
head(abf, 5) %>% 
  select(rank, snp, zscore, BF, ABF, ABFinv, P1_ABF) %>%
  pander
```


-------------------------------------------------------------------
 rank   snp    zscore      BF          ABF      ABFinv     P1_ABF  
------ ------ -------- ----------- ----------- --------- ----------
  1     rs23   5.773    8.692e-07   8.692e-07   1150477    0.7282  

  2     rs17   -5.585   2.566e-06   2.566e-06   389741     0.2467  

  3     rs42   -4.985   6.484e-05   6.484e-05    15423    0.009762 

  4     rs11   -4.846   8.525e-05   8.525e-05    11730    0.007425 

  5     rs26   4.718    0.0002211   0.0002211    4522     0.002862 
-------------------------------------------------------------------

# FINEMAP


```r
finemap <- finemapr(ss, ld, n = as.integer(1/0.015^2), 
  tool = "finemap", args = "--n-causal-max 1")
```

Top 5 SNPs:


```r
head(finemap$snp[[1]], 5) %>% 
  pander
```


--------------------------------------------------------------------
 snp    rank_z   rank_pp   snp_prob   snp_prob_cumsum   snp_log10bf 
------ -------- --------- ---------- ----------------- -------------
 rs23     1         1       0.7008        0.7009           2.042    

 rs17     2         2       0.2642        0.9651           1.227    

 rs42     3         3       0.0144        0.9795          -0.1637   

 rs11     4         4       0.0077        0.9872          -0.4397   

 rs26     5         5       0.0044        0.9916          -0.6838   
--------------------------------------------------------------------

# Final plot: ABF vs. FINEMAP


```r
# join results from two methods
tab <- full_join(
    select(finemap$snp[[1]], rank_z, snp, snp_prob) %>% rename(post_finemap = snp_prob),  
    select(abf, snp, P1_ABF) %>% rename(post_abf = P1_ABF),
    by = "snp") %>%
  head(9)

# table for plot
ptab <- gather(tab, method, post_prob, -rank_z, -snp) %>%
  mutate(
    snp_rank_z = paste0("#", rank_z, "\n", snp),
    method = factor(method, labels = c("ABF (Wakefield)", "FINEMAP (Benner)")))
  
# plot
ggplot(ptab, aes(factor(snp_rank_z), post_prob, fill = method)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "SNP", y = "Posterior probability")
```

![](figures/final_plot-1.png)<!-- -->
