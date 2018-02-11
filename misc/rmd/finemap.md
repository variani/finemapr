---
output:
  html_document:
    keep_md: true
---




```r
# load example data from FINEMAP (http://www.christianbenner.com/)
ex <- example_finemap()

# run the tool
# - simulated data in region 1 assumed 2 causal SNPs
out <- run_finemap(ex$tab1, ex$ld1, ex$n1, args = "--n-causal-max 3")
```


```r
print(out)
```

```
 - command: ~/apps/finemap/finemap --sss --log --n-causal-max 3 --in-files region.master 
 - config:
# A tibble: 7,506 x 4
   rank config         config_prob config_log10bf
  <int> <chr>                <dbl>          <dbl>
1     1 rs15,rs47           0.607            42.7
2     2 rs15,rs42,rs47      0.0326           43.1
3     3 rs15,rs34,rs47      0.0222           42.9
# ... with 7,503 more rows
```


```r
plot(out, label_size = 3, grid_ncol = 1)
```

![](figures/plot_finemap-1.png)<!-- -->

