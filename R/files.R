#-----------------
# File names
#-----------------

#' @rdname FinemaprFinemap
#' @export
filename_zscore.FinemaprFinemap <- function(x, locus = 1) paste0(locus, ".region.z") 

#' @rdname FinemaprFinemap
#' @export
filename_ld.FinemaprFinemap <- function(x, locus = 1) paste0(locus, ".region.ld")

#' @rdname FinemaprFinemap
#' @export
filename_snp.FinemaprFinemap <- function(x, locus) paste0(locus, ".region.snp")

#' @rdname FinemaprFinemap
#' @export
filename_config.FinemaprFinemap <- function(x, locus) paste0(locus, ".region.config")

#' @rdname FinemaprFinemap
#' @export
filename_log.FinemaprFinemap <- function(x, locus) paste0(locus, ".region.log")

#' @rdname FinemaprFinemap
#' @export
filename_master.FinemaprFinemap <- function(x) "region.master"
