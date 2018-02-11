#' S3 class FinemaprFinemap.
#'
#' @name FinemaprFinemap
#' @rdname FinemaprFinemap
#'
#' @exportClass FinemaprFinemap

#' @rdname FinemaprFinemap
#' @export
print.FinemaprFinemap <- function(x, ...)
{
  cat(" - command:", x$cmd, "\n")
  
  if(x$status) {
    cat(" - config:\n")
    print(x$config, n = 3)
  }
}

#' @rdname FinemaprFinemap
#' @export
plot.FinemaprFinemap <- function(x, 
  grid_nrow = NULL, grid_ncol = NULL, 
  label_size = getOption("finemapr_label_size"),
  label_size_config = label_size, label_size_snp = label_size,  
  top_rank = getOption("top_rank"),
  top_rank_config = top_rank, top_rank_snp = top_rank,    
  ...)
{
  p1 <- plot_ncausal(x, ...)
  p2 <- plot_config(x,  
    top_rank = top_rank_config, 
    label_size = label_size_config, ...)
  p3 <- plot_snp(x, 
    top_rank = top_rank_snp,
    label_size = label_size_snp, ...)
  
  plot_grid(p1, p2, p3, nrow = grid_nrow, ncol = grid_ncol)
}

#---------------
# Custom methods
#---------------

#' @export plot_ncausal
plot_ncausal <- function(x, ...) UseMethod("plot_ncausal")

#' @export plot_config
plot_config <- function(x, ...) UseMethod("plot_config")

#' @export plot_snp
plot_snp <- function(x, ...) UseMethod("plot_snp")

#' @rdname FinemaprFinemap
#' @export
plot_ncausal.FinemaprFinemap <- function(x, 
  ...)
{
  ptab <- x$ncausal
  
  sum_prop_zero <- filter(ptab, ncausal_num == 0)[["prob"]]  %>% sum
  if(sum_prop_zero == 0) {
    ptab <- filter(ptab, ncausal_num != 0)
  }
  
  ptab <- mutate(ptab, 
    ncausal_num = factor(ncausal_num, levels = sort(unique(ncausal_num), decreasing = TRUE)),
    type = factor(type, levels = c("prior", "post")))
    
  p <- ggplot(ptab, aes(ncausal_num, ncausal_prob, fill = type)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    coord_flip() + theme(legend.position = "top") + 
    scale_fill_manual(values = c("grey50", "orange"))
  
  return(p)
}

#' @rdname FinemaprFinemap
#' @export
plot_config.FinemaprFinemap <- function(x, xlim_config = c(0, 1.5), 
  label_size = getOption("finemapr_label_size"),  
  top_rank = getOption("top_rank"),  
  ...)
{
  ptab <- x$config

  ptab <- head(ptab, top_rank)

  ptab <- mutate(ptab,
    label = paste0(config, "\n", 
      "P = ", round(config_prob, 2),
      "; ", "log10(BF) = ", round(config_log10bf, 2)))

  ggplot(ptab, aes(config_prob, rank)) + 
    geom_point() + 
    geom_segment(aes(xend = config_prob, yend = rank, x = 0)) + 
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.025, size = label_size) + 
    xlim(xlim_config) + 
    scale_y_continuous(limits  = c(top_rank + 0.5, 0.5), trans = "reverse")
}

#' @rdname FinemaprFinemap
#' @export
plot_snp.FinemaprFinemap <- function(x, xlim_snp = c(0, 1.5), 
  label_size = getOption("finemapr_label_size"),  
  top_rank = getOption("top_rank"),  
  ...)
{
  ptab <- x$snp

  ptab <- head(ptab, top_rank)

  ptab <- mutate(ptab,
    rank = seq(1, n()), 
    label = paste0(snp, "\n", 
      "P = ", round(snp_prob, 2),
      "; ", "log10(BF) = ", round(snp_log10bf, 2)))

  ggplot(ptab, aes(snp_prob, rank)) + 
    geom_point() + 
    geom_segment(aes(xend = snp_prob, yend = rank, x = 0)) + 
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.025, size = label_size) + 
    xlim(xlim_snp) + 
    scale_y_continuous(limits  = c(top_rank + 0.5, 0.5), trans = "reverse")
}
