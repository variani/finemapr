
#' @export
ref_ld <- function(chr = 16, start = 53.767e6, end = 53.768e6, mono = FALSE,
  pop = "EUR", snps = NULL,
  measure = c("r2", "r", "D"), tol = 1e-10,
  strict = FALSE,
  verbose = 0)
{
  measure <- match.arg(measure)
  
  stopifnot(requireNamespace("gaston"))
    
  bed <- ref_bed(chr = chr, start = start, end = end, mono = mono,
  pop = pop, snps = snps,
  strict = strict, verbose = verbose)

  ld <- gaston::LD(bed, c(1, ncol(bed)), measure = measure)
  
  ld[abs(ld) < tol] <- 0
  
  return(ld)
}

#' @export
ref_pop <- function()
{
  path.package("finemapr") %>% 
    file.path("inst/extdata/pop_label_phase3_2504_without_rels") %>%
    read_delim(delim = " ", col_types = cols())
}

#' @export
ref_bed <- function(chr = 16, start = 53.767e6, end = 53.768e6, mono = FALSE,
  pop = "EUR", snps = NULL,
  strict = FALSE,
  verbose = 0)
{
  ### inc
  stopifnot(requireNamespace("RCurl"))
  stopifnot(requireNamespace("gaston"))
  
  ### agrs
  groups_pop <- pop # rename  
  
  ### url
  url <- paste0("http://tabix.iobio.io/?cmd=-h%20%27",
    "http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/",
    "vcf.b37/",
    "chr", chr, ".1kg.phase3.v5a.vcf.gz",
    "%27%20", chr, ":", start, "-", end)
  
  if(verbose) {
    cat(" - url:", url, "\n")
  }
  
  ### download vcf
  vcf <- RCurl::getURL(url)
  
  ### convert from vcf to bed (plink)
  file_tmp <- tempfile()
  write_file(vcf, file_tmp)
  
  bed <- gaston::read.vcf(file_tmp)
  ret <- unlink(file_tmp)
  
  ### filter by population labels
  pop <- ref_pop() 
  ids <- with(pop, sample[super_pop %in% groups_pop])
  stopifnot(all(ids %in% bed@ped$id))
  
  if(verbose) {
    cat(" - #individuals:", nrow(bed), "\n")
  }
  bed <- bed[bed@ped$id %in% ids, ]
  if(verbose) {
    cat("  -- #individuals after pop. filter:", nrow(bed), "\n")
  }

  ### filter monomorphic snps
  if(verbose) {
    cat(" - #snps:", ncol(bed), "\n")
  }
  if(!mono) {
    bed <- bed[, bed@sigma > 0]
    if(verbose) {
      cat("  -- #snps after mono. filter:", ncol(bed), "\n")
    }
  }
  
  
  ### filter by `snps`
  if(!is.null(snps)) {
    stopifnot(any(snps %in% bed@snps$id))
    
    if(!all(snps %in% bed@snps$id)) {
      msg <- paste("not all `snps` are found")
      
      if(strict) {
        stop(msg)
      } else {
        warning(msg)
      }
    }
    
    bed <- bed[, bed@snps$id %in% snps]
  }

  return(bed)
}

#----------------------------
# Functions from proxysnps
#----------------------------

#'
#' https://github.com/slowkow/proxysnps/blob/master/R/get_vcf.R
#'
#' @export
proxysnps_ref_vcf <- function(chr = 16, start = 53.767e6, end = 53.768e6, pop = "EUR",
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
proxysnps_ref_ld <- function(vcf, tol = 1e-10, measure = c("R.squared", "D.prime"),
  verbose = 0)
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

