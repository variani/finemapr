#----------------------
# Class Declatation
#----------------------

#' S3 class Finemapr.
#'
#' @name Finemapr
#' @rdname Finemapr
#'
#' @exportClass Finemapr


#----------------------
# Methods Declatation
#----------------------

#' @export plot_ncausal
plot_ncausal <- function(x, ...) UseMethod("plot_ncausal")

#' @export plot_config
plot_config <- function(x, ...) UseMethod("plot_config")

#' @export plot_snp
plot_snp <- function(x, ...) UseMethod("plot_snp")

#' @export plot_zscore
plot_zscore <- function(x, ...) UseMethod("plot_zscore")

#' @export process_tab
process_tab <- function(x, ...) UseMethod("process_tab")

#' @export process_ld
process_ld <- function(x, ...) UseMethod("process_ld")

#' @export process_n
process_n <- function(x, ...) UseMethod("process_n")

#' @export process_annot
process_annot <- function(x, ...) UseMethod("process_annot")

#' @export write_files
write_files <- function(x, ...) UseMethod("write_files")

#' @export run_tool
run_tool <- function(x, ...) UseMethod("run_tool")

#' @export collect_results
collect_results <- function(x, ...) UseMethod("collect_results")

#' @export extract_credible_set
extract_credible_set <- function(x, ...) UseMethod("extract_credible_set")

#------------------------
# File names
#------------------------

#' @export filename_zscore
filename_zscore <- function(x, ...) UseMethod("filename_zscore")

#' @export filename_ld
filename_ld <- function(x, ...) UseMethod("filename_ld")

#' @export filename_snp
filename_snp <- function(x, ...) UseMethod("filename_snp")

#' @export filename_config
filename_config <- function(x, ...) UseMethod("filename_config")

#' @export filename_k
filename_k <- function(x, ...) UseMethod("filename_k")

#' @export filename_log
filename_log <- function(x, ...) UseMethod("filename_log")

#' @export filename_master
filename_master <- function(x, ...) UseMethod("filename_master")

#' @export filename_annot
filename_annot <- function(x, ...) UseMethod("filename_annot")

#------------------------
# Print methods
#------------------------

#' @rdname Finemapr
#' @export
print.Finemapr <- function(x, ...)
{
  cat(" - command:", x$cmd, "\n")
    
  if(x$status) {
    cat(" - snp:\n")
    print(x$snp, n = 3)
  }
}


#------------------------
# Plot methods
#------------------------

#' @rdname Finemapr
#' @export
print.Finemapr <- function(x, ...)
{
  cat(" - cmd:", x$cmd, "\n")
  cat(" - #loci:", x$num_loci, "\n") 
  
  ret <- lapply(seq(1, x$num_loci), function(i) {    
    cat(" - locus:",i, "\n")
    cat("  -- snp:\n")
    print(x$snp[[i]], n = 3)
    cat("  -- ", length(x$snps_credible[[i]]), " snps in ",
      100*x$prop_credible, "% credible set", 
      ": ", paste(x$snps_credible[[i]], collapse = ", "), "...", 
      "\n", sep = "") 
  })

  return(invisible())
}  
  
#' @rdname Finemapr
#' @export
plot.Finemapr <- function(x, ...)
{
  plot_snp(x, ...)
}


#' @rdname Finemapr
#' @export
plot_snp.Finemapr <- function(x, lim_prob = c(0, 1.5), 
  label_size = getOption("finemapr_label_size"),  
  top_rank = getOption("top_rank"),  
  ...)
{
  ptab <- x$snp

  ptab <- head(ptab, top_rank)

  ptab <- mutate(ptab,
    label = paste0(snp, "\n", 
      "P = ", round(snp_prob, 2)))
      
  ggplot(ptab, aes(snp_prob, rank)) +
    geom_vline(xintercept = 1, linetype = 3) + 
    geom_point() + 
    geom_segment(aes(xend = snp_prob, yend = rank, x = 0)) + 
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.025, size = label_size) + 
    xlim(lim_prob) + 
    scale_y_continuous(limits  = c(top_rank + 0.5, 0.5), trans = "reverse")
}

#------------------------
# Other methods
#------------------------

#' @rdname Finemapr
#' @export
extract_credible_set.Finemapr <- function(x, ...)
{
  lapply(x$snp, function(snp) {
    snp_below <- snp %>% filter(snp_prob_cumsum <= x$prop_credible)
    snps <- head(snp, nrow(snp_below) + 1) %$% snp
    
    # the case: the single top snps covers 100% of credibility
    #if(length(snps) == 0) {
    #  snps <- head(snp, 1) %$% snp
    #}
    
    return(snps)
  })
}



