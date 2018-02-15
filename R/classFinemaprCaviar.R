#' S3 class FinemaprCaviar.
#'
#' @name FinemaprCaviar
#' @rdname FinemaprCaviar
#'
#' @exportClass FinemaprCaviar

#' @rdname FinemaprCaviar
#' @export
print.FinemaprCaviar <- function(x, ...)
{
  cat(" - command:", x$cmd, "\n")
    
  if(x$status) {
    #cat(" - see log output in `log`\n")
    cat(" - tables of results: `snp`\n")
    
    cat(" - snp:\n")
    print(x$snp, n = 3)

    cat(" - 95%-causal set (ordered):", paste(x$set, collapse = ", "), "\n")
  }
}

#' @rdname FinemaprCaviar
#' @export
plot.FinemaprCaviar <- function(x, ...)
{
  plot_snp(x, ...)
}

#---------------
# Custom methods
#---------------

#' @rdname FinemaprCaviar
#' @export
plot_snp.FinemaprCaviar <- function(x, lim_prob = c(0, 1.5), 
  label_size = getOption("finemapr_label_size"),  
  top_rank = getOption("top_rank"),  
  ...)
{
  ptab <- x$snp

  ptab <- head(ptab, top_rank)

  ptab <- mutate(ptab,
    label = paste0(snp, "\n", 
      "P = ", round(snp_prob, 2),
      "; ", "P(set) = ", round(snp_prob_set, 2)))

  ggplot(ptab, aes(snp_prob, rank)) +
    geom_vline(xintercept = 1, linetype = 3) + 
    geom_point() + 
    geom_segment(aes(xend = snp_prob, yend = rank, x = 0)) + 
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.025, size = label_size) + 
    xlim(lim_prob) + 
    scale_y_continuous(limits  = c(top_rank + 0.5, 0.5), trans = "reverse")
}
