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
    cat(" - config (top 5):\n")
    print(head(x$config, 5))
  }
}

#' @rdname FinemaprFinemap
#' @export
plot.FinemaprFinemap <- function(x, ...)
{
  ptab <- x$causal
  
  sum_prop_zero <- filter(tab, num == 0)[["prob"]]  %>% sum
  if(sum_prop_zero == 0) {
    ptab <- filter(ptab, num != 0)
  }
  
  ptab <- mutate(ptab, 
    num = factor(num, levels = sort(unique(num), decreasing = TRUE)),
    type = factor(type, levels = c("prior", "post")))
    
  p <- ggplot(ptab, aes(num, prob, fill = type)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    coord_flip() + theme(legend.position = "top") + 
    scale_fill_manual(values = c("grey50", "orange"))
  
  return(p)
}

