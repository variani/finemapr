finemapr_guess_finemap <- function()
{
  switch(Sys.info()[["nodename"]],
    "topoki" = "~/apps/finemap/finemap",
    "finemap")
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  
  op_finemapr <- list(
    finemapr_finemap = finemapr_guess_finemap()
  )
  
  ind_set <- !(names(op_finemapr) %in% names(op))
  if(any(ind_set)) {
    options(op_finemapr[ind_set])
  }
  
  invisible()
}
