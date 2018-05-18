#------------------------
# Plot methods
#------------------------


#' @rdname Finemapr
#'
#' @note
#' Colors: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
#'
#' @export
plot_zscore.Finemapr <- function(x, locus = 1, 
  label_size = getOption("finemapr_label_size"), 
  selected,
  # ggrepel
  force = 20,
  color_main = "grey40", color_selected = "aquamarine4", color_credible = "dodgerblue4", color_missing = "brown",
  ...)
{
  ### arg
  missing_selected <- missing(selected)
  
  ### data
  tab <- x$tab[[locus]]
  tab <- mutate(tab,
    pval = pchisq(zscore^2, df = 1, lower.tail = FALSE))
      
  p <- ggplot(tab, aes(pos, -log10(pval))) + geom_point(color = color_main)
  
  ### sel snps
  if(!missing_selected) {
    tab_selected <- filter(tab, snp %in% selected)
    #(nrow(tab_selected) == length(selected))
  
    if(nrow(tab_selected)) {
      p <- p + 
        geom_point(
          data = tab_selected, aes(pos, -log10(pval)), color = color_selected) +
        geom_text_repel(
          data = tab_selected, aes(label = snp),
          force = force, size = label_size, color = color_selected)
    } 
  }
  
  ### credible set
  if(!is.null(x$snps_credible)) {
    credible <- x$snps_credible[[locus]]
    
    tab_credible <- filter(tab, snp %in% credible)
    stopifnot(nrow(tab_credible) == length(credible))
  
    p <- p + 
      geom_point(
        data = tab_credible, aes(pos, -log10(pval)), color = color_credible) 
    
    tab_credible_top <- filter(tab, snp %in% head(credible, 3))
    p <- p + 
      geom_text_repel(
        data = tab_credible_top, aes(label = snp),
        force = force, size = label_size, color = color_credible)
  }

  ### missing snps 
  if(!is.null(x$snps_missing_finemap)) {
    snps_missing <- x$snps_missing_finemap[[locus]]
    
    tab_missing <- filter(tab, snp %in% snps_missing)
    stopifnot(nrow(tab_missing) == length(snps_missing))
  
    p <- p + 
      geom_point(
        data = tab_missing, aes(pos, -log10(pval)), color = color_missing) 
    
    tab_missing_top <- head(tab_missing, 3) 
    p <- p + 
      geom_text_repel(
        data = tab_missing_top, aes(label = snp),
        force = force, size = label_size, color = color_missing)
  }
  
  return(p)    
}
