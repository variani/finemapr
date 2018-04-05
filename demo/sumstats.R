### inc
library(dplyr)

library(rsnps)

### get data on summary stats from GWAS of BMI
if(!exists("ss")) {
  ss <- sumstats_bmi()
}

ss %>% arrange(P) %>% head(5)

ncbi_snp_query("rs1421085") # chr16:53,800,954 

### select SNPs in a region of the top SNP rs1421085
# - not possible, the data of `ss` has no info. about chr/pos
# ss_region <- filter(ss, ...)

### get ref. genotypes
if(!exists("vcf")) {
  vcf <- ref_vcf(chr = 16, start = 53.800e6, end = 53.801e6, pop = "EUR")
}

### compute LD matrix
ld <- ref_ld(vcf)

image(ld)

### overlap snps
ss_region <- data_frame(SNP = rownames(ld)) %>% 
  left_join(ss, by = "SNP") %>%
  filter(!is.na(Z))
  
snps_region <- ss_region$SNP

ld_region <- ld[snps_region, snps_region]

### run fine-mapping analysis 
zscore_region <- ss_region %>% select(SNP, Z)

out <- run_paintor(zscore_region, as.matrix(ld_region), n_avr, args = "-enumerate 3")

plot(out)

### additional example: get ref. genotypes for a given set of snps
if(!exists("vcf_snps")) {
  vcf_snps <- ref_vcf(chr = 16, start = 53.800e6, end = 53.801e6, pop = "EUR",
    snps = c("rs12446228", "rs9939973", "rs9940128", "rs1421085"))
    
  # the computation of LD is faster on a reduced set of snps
  ld_snps <- ref_ld(vcf_snps)
}


