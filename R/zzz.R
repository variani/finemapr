finemapr_guess_finemap <- function()
{
  switch(Sys.info()[["nodename"]],
    "topoki" = "~/apps/finemap/finemap",
    "finemap")
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  
  op_finemapr <- list(
    # paths to tools
    finemapr_finemap = finemapr_guess_finemap(),
    # plot
    finemapr_label_size = 4,
    top_rank = 5
  )
  
  ind_set <- !(names(op_finemapr) %in% names(op))
  if(any(ind_set)) {
    options(op_finemapr[ind_set])
  }
  
  invisible()
}
