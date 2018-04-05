#'
#' https://github.com/slowkow/proxysnps/blob/master/R/get_vcf.R
#'
#' @export
ref_vcf <- function(chr = 16, start = 53.767e6, end = 53.768e6, pop = "EUR",
  snps = NULL,
  verbose = 0)
{
  ### inc
  stopifnot(requireNamespace("proxysnps"))
  
  ### get genotypes from the ref. panel
  if(verbose) {
    cat(" - get reference genotypes\n")
  }  
  vcf <- proxysnps::get_vcf(chrom = chr, start = start, end = end, pop = pop)

  ### filter
  if(!is.null(snps)) {
    vcf$meta <- filter(vcf$meta, ID %in% snps)
    
    names_snp <- rownames(vcf$geno)
    ind_snps <- names_snp %in% snps
    vcf$geno <- vcf$geno[ind_snps, ]
  }
  
  return(vcf)
}

#'
#' https://github.com/slowkow/proxysnps/blob/master/R/get_vcf.R
#'
#' @export
ref_ld <- function(vcf, tol = 1e-10, measure = c("R.squared", "D.prime"),
  verbose = 1)
{
  ### args
  stopifnot(!missing(vcf))
  measure <- match.arg(measure)
  
  ### inc
  stopifnot(requireNamespace("proxysnps"))
  
  ### the ref. panel
  geno <- vcf$geno
  
  num_snps <- nrow(geno)
  names_snps <- rownames(geno)
  
  ### build LD matrix
  ld_elem <- lapply(seq(1, num_snps), function(i) {
    if(verbose) {
      cat(" - compute ld for snp:", i, "/", num_snps, "\n")
    }
    
    j_val <- seq(i, num_snps)
    list(i = rep(i, length(j_val)), j = j_val, 
      x = sapply(j_val, function(j) {
        ifelse(i == j, 1, proxysnps::compute_ld(geno[i, ], geno[j, ])[[measure]])
      }))
  })
  
  ld <- sparseMatrix(
    i = sapply(ld_elem, function(x) x$i) %>% unlist,
    j = sapply(ld_elem, function(x) x$j) %>% unlist,
    x = sapply(ld_elem, function(x) x$x) %>% unlist,
    symmetric = TRUE,
  )
  
  ld <- drop0(ld, tol = tol) 

  rownames(ld) <- names_snps
  colnames(ld) <- names_snps

  return(ld)  
}

