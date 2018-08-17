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
  color_main = "grey75", color_selected = "aquamarine4", 
  color_credible = "dodgerblue4", color_missing = "brown",
  plot_missing = FALSE,
  ...)
{
  ### arg
  missing_selected <- missing(selected)
  
  ### data
  tab <- x$tab[[locus]]
  tab <- mutate(tab,
    pval = pchisq(zscore^2, df = 1, lower.tail = FALSE))
  
  if(!plot_missing) {
    tab <- filter(tab, finemap)
  }
  
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
    
    tab_credible_top <- filter(tab, snp %in% head(credible, 10))
    p <- p + 
      geom_text_repel(
        data = tab_credible_top, aes(label = snp),
        force = force, size = label_size, color = color_credible)
  }

  ### missing snps 
  if(!is.null(x$snps_missing_finemap) & plot_missing) {
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
  
  ### labs
  p <- p + 
    scale_x_continuous(labels = scales::comma) +
    labs(x = "Position (bp)", y = expression(-log[10](P)))
  
  ### cleaner theme
  p <- p + theme(panel.grid.minor = element_blank())
  
  return(p)    
}
