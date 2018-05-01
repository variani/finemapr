
#' Run FINEMAP.
#'
#' @examples
#' ex <- example_finemap()
#' out <- finemapr(ex$tab1, ex$ld1, ex$n1, args = "--n-causal-max 1")
#'
#' @export
finemapr <- function(tab, ld, n, 
  method = c("finemap"),
  dir_run,
  tool = getOption("finemapr_finemap"), args = "")
{
  ### arg
  method <- match.arg(method)
  
  stopifnot(!missing(tab))
  stopifnot(!missing(ld))
  stopifnot(!missing(n))
  
  ### create an object of class `Finemapr`: basic slots and class attribute
  out <- list(method = method, 
    dir_run = paste("run", method, sep = "_"), args = args,
    num_loci = ifelse(class(tab)[1] == "list", length(tab), 1))

  class_finemapr <- switch(out$method,
    "finemap" = "FinemaprFinemap",
    stop("switch error on `method`"))
  oldClass(out) <- c(class_finemapr, "Finemapr", oldClass(out))

  ### process input
  out <- process_tab(out, tab)
  out <- process_ld(out, ld)
  # out <- process_n(out, n)
  
  ### write files 
  #write_files(out)
  
  #### run
  #ret <- run_tool(out)

  #### read results
  #out <- read_results(out)
    
  return(out)
}

#' @rdname Finemapr
#' @export
process_tab.Finemapr <- function(x, tabs, ...)
{
  ### process input
  if(class(tabs)[1] != "list") {
    tabs <- list(tabs)
  }
  
  stopifnot(length(tabs) == x$num_loci)
  
  ### prepare tabels of Z-scores
  out_tabs <- lapply(tabs, function(tab) {
    tab <- as_data_frame(tab)

    stopifnot(ncol(tab) >= 2)
    
    names_all <- names(tab)
    names_select <- c(
      finemapr_find_name("snp", names_all, strict = TRUE),
      finemapr_find_name("zscore", names_all, strict = TRUE))
    
    tab <- select_(tab, .dots = names_select)
    names(tab) <- finemapr_names_tab()
    
    # manage missing Z-scores
    snps_zscore_missing <- filter(tab, is.na(zscore)) %$% snp 
    tab <- filter(tab, !is.na(zscore)) 
  
    list(tab = tab, 
      snps_zscore_missing = snps_zscore_missing)
  })
  
  ### write back to `x` and return
  x$tab <- lapply(out_tabs, function(x) x$tab)
  x$snps_zscore <- lapply(out_tabs, function(x) x$tab[[finemapr_names_tab_snp()]])
    
  x$snps_zscore_missing <- lapply(out_tabs, function(x) x$snps_zscore_missing)
  
  return(x)
}

#' @rdname Finemapr
#' @export
process_ld.Finemapr <- function(x, lds, ...)
{
  ### process input
  if(class(lds)[1] != "list") {
    lds <- list(lds)
  }
  
  stopifnot(length(lds) == x$num_loci)
  
  ### prepare tabels of Z-scores
  out_lds <- lapply(seq_along(lds), function(locus) {
    ld <- lds[[locus]]
    
    stopifnot(class(ld) == "matrix")
    stopifnot(!is.null(colnames(ld)))
    stopifnot(!is.null(rownames(ld)))
    
    # manage SNP names across variables: ld, zscore    
    snps_ld <- colnames(ld)
    snps_zscore <- x$snps_zscore[[locus]]
    
    ind <- snps_ld %in% snps_zscore
    snps_finemap <- snps_ld[ind]
    snps_ld_missing <- snps_ld[!ind]
    
    # check the proportion of `snps_ld_missing`
    prop_snps_missing <- length(snps_ld_missing) / 
      (length(snps_ld_missing) + length(snps_finemap))
    stopifnot(prop_snps_missing < 0.20)
    
    # subset LD matrix
    ld <- ld[snps_finemap, snps_finemap]
    
    # some tools require all diagonals to be `1`
    stopifnot(all(round(diag(ld), 4) == 1))
    diag(ld) <- 1
   
    list(ld = ld, 
      snps_ld_missing = snps_ld_missing)
  })
  
  x$ld <- lapply(out_lds, function(x) x$ld)
  x$snps_ld_missing <- lapply(out_lds, function(x) x$snps_ld_missing)  
  x$snps_finemap <- lapply(out_lds, function(x) colnames(x$ld))
  
  return(x)
}


