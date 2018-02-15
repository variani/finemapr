






```r
# specify path to PAINTOR (user-specific)
options(finemapr_paintor = "~/apps/paintor/PAINTOR")

# load the same example data as for FINEMAP
ex <- example_finemap()

# run the tool
# - arg. #1: table with SNP name and Z-score in the first 2 columns
# - arg. #2: the LD matrix with colnames/rownames corresponding to SNP names
# - arg. #3 (optional; the default value is 1,000,000): the number of invidivudals
# - arg. `args`: pass other arguments to the tool as a string 
out_paintor <- run_paintor(ex$tab1, ex$ld1, ex$n1, args = "-enumerate 3")
```


```r
print(out_paintor)
```

```
 - command: ~/apps/paintor/PAINTOR -input region.master -Zhead zscore -LDname ld -annotations dummy_ones -enumerate 3 -num_samples 5363 
 - tables of results: `snp`
 - snp:
# A tibble: 50 x 4
   rank   snp     zscore snp_prob
  <int> <chr>      <dbl>    <dbl>
1     1  rs15 10.9151878 1.000000
2     2  rs47  8.7995877 1.000000
3     3  rs25 -0.1326266 0.020042
# ... with 47 more rows
 - annotations: dummy_ones 
 - logBF (proportional to the model likelihood): 100.7274 
```


```r
plot(out_paintor, label_size = 3)
```

![](figures/plot_paintor-1.png)<!-- -->

