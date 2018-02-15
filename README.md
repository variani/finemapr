# finemapr

R interface to fine-mappers:

- FINEMAP http://www.christianbenner.com/
- CAVIAR https://github.com/fhormoz/caviar
- PAINTOR https://github.com/gkichaev/PAINTOR_V3.0

By using `finemapr`, your input files are automatically prepared for each tool, the analysis workflow is tool-independent; and exploration of fine-mapping results is powered by R in printing/plotting/data export.

## Unified analysis workflow

```r
# set up
options(finemapr_finemap = "<path to fine-mapping tool>")

# read input files
my_zscores <- read_zscores("<my_scores.tab>")
my_ld <- read_ld("<my_ld.tab>")

# run analysis
out <- run_<tool>(my_zscores, my_ld, args = "<custom arguments>")

# explore results
print(out)
head(out$snp) # main table of results
plot(out)

# export results
write.table(out$snp, "<my_results.tab>")
```

## FINEMAP example

```r
# specify where the FINEMAP tool is located (user-specific)
# - the user needs to download and install the tool before the analysis
# - here, the tool binary file is stored in `~/apps/finemap/` directory
options(finemapr_finemap = "~/apps/finemap/finemap")

# load example data from FINEMAP (http://www.christianbenner.com/)
# - data are obtained from the same directory `~/apps/finemap/`, subdirectory `example`
ex <- example_finemap()

# run the tool
# - simulated data in region 1 (used here in the analysis) assumed 2 causal SNPs
# - arg. #1: table with SNP name and Z-score in the first 2 columns
# - arg. #2: the LD matrix with colnames/rownames corresponding to SNP names
# - arg. #3: the number of invidivudals
# - arg. `args`: pass other arguments to the tool as a string 
out <- run_finemap(ex$tab1, ex$ld1, ex$n1, args = "--n-causal-max 3")
```


```r
print(out)
```

```
 - command: ~/apps/finemap/finemap --sss --log --n-causal-max 3 --in-files region.master 
 - see log output in `log`
 - tables of results: `config`, `snp`, `ncausal`
 - config:
# A tibble: 8,410 x 4
   rank config         config_prob config_log10bf
  <int> <chr>                <dbl>          <dbl>
1     1 rs15,rs47           0.607            42.7
2     2 rs15,rs42,rs47      0.0326           43.1
3     3 rs15,rs34,rs47      0.0222           42.9
# ... with 8,407 more rows
```


```r
plot(out, label_size = 3, grid_ncol = 1)
```

![](misc/rmd/figures/plot_finemap-1.png)


## CAVIAR example

```r
# specify path to CAVIAR (user-specific)
options(finemapr_caviar = "~/apps/caviar/CAVIAR")

# load the same example data as for FINEMAP
ex <- example_finemap()

# run the tool
# - arg. #1: table with SNP name and Z-score in the first 2 columns
# - arg. #2: the LD matrix with colnames/rownames corresponding to SNP names
# - arg. `args`: pass other arguments to the tool as a string 
out_caviar <- run_caviar(ex$tab1, ex$ld1, args = "-c 3")
```


```r
print(out_caviar)
```

```
 - command: ~/apps/caviar/CAVIAR -c 3 -z region.z -l region.ld -o log 
 - tables of results: `snp`
 - snp:
# A tibble: 50 x 4
   rank snp   snp_prob_set snp_prob
  <int> <chr>        <dbl>    <dbl>
1     1 rs15        0.439    1.00  
2     2 rs47        0.439    1.00  
3     3 rs42        0.0120   0.0274
# ... with 47 more rows
 - 95%-causal set (ordered): rs15, rs47, rs42, rs34, rs20, rs5, rs27, rs45, rs38, rs17, rs25, rs18, rs11, rs19, rs44, rs40, rs8, rs24 
```


```r
plot(out_caviar, label_size = 3)
```

![](misc/rmd/figures/plot_caviar-1.png)
