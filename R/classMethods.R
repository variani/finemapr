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

#' @export process_tab
process_tab <- function(x, ...) UseMethod("process_tab")

#' @export process_ld
process_ld <- function(x, ...) UseMethod("process_ld")

#' @export process_n
process_n <- function(x, ...) UseMethod("process_n")

#' @export write_files
write_files <- function(x, ...) UseMethod("write_files")

#' @export run_tool
run_tool <- function(x, ...) UseMethod("run_tool")

#' @export collect_results
collect_results <- function(x, ...) UseMethod("collect_results")

#------------------------
# File names
#------------------------

#' @export filename_
filename_zscore <- function(x, ...) UseMethod("filename_zscore")

#' @export filename_
filename_ld <- function(x, ...) UseMethod("filename_ld")

#' @export filename_
filename_snp <- function(x, ...) UseMethod("filename_snp")

#' @export filename_
filename_config <- function(x, ...) UseMethod("filename_config")

#' @export filename_
filename_log <- function(x, ...) UseMethod("filename_log")

#' @export filename_
filename_master <- function(x, ...) UseMethod("filename_master")


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
