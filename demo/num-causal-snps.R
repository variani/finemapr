# demo: shows how to specify the number of causal loci
# - finemap: note that the argument name `--n-causal-max`/`--n-causal-snps` 
#   differs by finemap version (http://www.christianbenner.com/)

### load example data
ex <- example_finemap()

### run finemap
# `--n-causal-max` is missing; the number of causal snps is infered from `prior_k`
out1 <- finemapr(ex$tab1, ex$ld1, ex$n1, method = "finemap", prior_k = rep(1, 4))

# values specifed in `--n-causal-max` and `prior_k` are concordant
# (the number of causal snps = 2)
out2 <- finemapr(ex$tab1, ex$ld1, ex$n1, method = "finemap", prior_k = rep(1, 2), args = "--n-causal-max 2")

# values specifed in `--n-causal-max` and `prior_k` are in conflict, but `prior_k` matters
# (the number of causal snps = 2)
out3 <- finemapr(ex$tab1, ex$ld1, ex$n1, method = "finemap", prior_k = rep(1, 2), args = "--n-causal-max 4")

# (the number of causal snps = 3) specified in `--n-causal-max`
out4 <- finemapr(ex$tab1, ex$ld1, ex$n1, method = "finemap",  args = "--n-causal-max 3")

