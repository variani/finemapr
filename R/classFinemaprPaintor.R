#' S3 class FinemaprPaintor.
#'
#' @name FinemaprPaintor
#' @rdname FinemaprPaintor
#'
#' @exportClass FinemaprPaintor

#' @rdname FinemaprPaintor
#' @export
print.FinemaprPaintor <- function(x, ...)
{
  cat(" - command:", x$cmd, "\n")
    
  if(x$status) {
    #cat(" - see log output in `log`\n")
    cat(" - tables of results: `snp`\n")
    
    cat(" - snp:\n")
    print(x$snp, n = 3)
  
    cat(" - annotations:", paste(x$annotations, collapse = ", "), "\n")
    cat(" - logBF (proportional to the model likelihood):", x$logBF, "\n")
  }
}


