#' Read z-scores from a file.
#'
#' @examples
#' f1 <- system.file("extdata/region1.z", package = "finemapr")
#' ztab <- read_zscore(f1)
#'
#' @export
read_zscore <- function(file_zscore)
{
  ### args
  stopifnot(!missing(file_zscore))
  stopifnot(file.exists(file_zscore))
  
  ### inc
  stopifnot(requireNamespace("data.table"))

  ### read
  tab <- data.table::fread(file_zscore) %>% as_data_frame
  
  stopifnot(ncol(tab) >= 2)
  
  ### processing column names
  # scenario 1: first two columns are `snp` and `zscore`
  tab <- tab[, 1:2]
  
  col1_class <- class(tab[[1]])
  col2_class <- class(tab[[2]])
  
  stopifnot(col1_class == "character")
  stopifnot(col2_class == "numeric")
  
  names(tab) <- c("snp", "zscore")
    
  return(tab)
}

#' Read LD matrix from a file.
#'
#' @examples
#' f1 <- system.file("extdata/region1.ld", package = "finemapr")
#' ld <- read_ld(f1)
#'
#' @export
read_ld <- function(file_ld, snps)
{
  ### args
  stopifnot(!missing(file_ld))
  stopifnot(file.exists(file_ld))
  
  ### vars
  missing_snps <- missing(snps)
  
  ### inc
  stopifnot(requireNamespace("data.table"))

  ### read
  mat <- data.table::fread(file_ld) %>% as.matrix
  
  stopifnot(ncol(mat) == nrow(mat))
  
  if(!missing_snps) {
    stopifnot(nrow(mat) == length(snps))
    
    rownames(mat) <- snps
    colnames(mat) <- snps
  }    
  
  return(mat)
}
