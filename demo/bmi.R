### par
window <- 500e3

### data
ss <- sumstats_bmi()

### region
# rs1421085: Chr16: 53800954 (on Assembly GRCh37)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5025863/

ld <- ref_ld(chr = 16, start = 53800954 - window/2, end = 53800954 + window/2, verbose = 1)

snps_region <- rownames(ld)
ss_region <- filter(ss, SNP %in% snps_region)

### fine mapping
n <- with(ss_region, round(mean(N), 0))

out <- finemapr(ss_region, ld, n, args = "--n-causal-max 3")

print(out)

plot(out)
