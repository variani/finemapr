library(devtools)
load_all("~/git/variani/finemapr")

library(magrittr)
library(dplyr)  

library(ggplot2)  
theme_set(theme_linedraw())

### data
file1_z <- system.file("extdata/region1.z", package = "finemapr")
file1_ld <- system.file("extdata/region1.ld", package = "finemapr")

z1 <- read_zscore(file1_z)
ld1 <- read_ld(file1_ld, snps = z1$snp)
n1 <- 5363

### run finemap
options(finemapr_finemap = "~/apps/finemap/finemap")

out_finemap <- run_finemap(z1, ld1, n1, args = "--n-causal-max 3")

### plot
plot(out_finemap, label_size = 3, 
  lim_prob_config = c(0, 1.6), lim_prob_snp = c(0, 3), lim_prob_ncausal = c(0, 1),
  grid_nrow = 1)

### save fig.
ggsave("misc/figures/finemap-2causal.png",
  width = 9, height = 4, units = "in", dpi = 200)
       
