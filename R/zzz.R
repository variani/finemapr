finemapr_guess_finemap <- function()
{
  switch(Sys.info()[["nodename"]],
    "topoki" = "~/apps/finemap/finemap",
    "tau.local" = "~/apps/finemap/finemap",
    "finemap")
}

finemapr_guess_caviar <- function()
{
  switch(Sys.info()[["nodename"]],
    "topoki" = "~/apps/caviar/CAVIAR",
    "tau.local" = "~/apps/caviar/CAVIAR",
    "finemap")
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  
  op_finemapr <- list(
    # paths to tools
    finemapr_finemap = finemapr_guess_finemap(),
    finemapr_caviar = finemapr_guess_caviar(),
    # plot
    finemapr_label_size = 4,
    top_rank = 5,
    lim_prob = c(0, 1.5)
  )
  
  ind_set <- !(names(op_finemapr) %in% names(op))
  if(any(ind_set)) {
    options(op_finemapr[ind_set])
  }
  
  invisible()
}
