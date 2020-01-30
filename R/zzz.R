finemapr_guess_tool <- function(tool = c("finemap", "caviar", "paintor", "cojo"))
{
  tool <- match.arg(tool)
  
  tool_bin <- switch(tool,
    "finemap" = "finemap",
    "caviar" = "CAVIAR",
    "paintor" = "paintor",
    "cojo" = "gcta",
    stop())
  
  path_apps <- file.path("~/.local/apps/", tool)
  file_apps <- file.path(path_apps, tool_bin)

  if(file.exists(file_apps)) {
    return(file_apps)
  }
  
  return(tool_bin)
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  
  op_finemapr <- list(
    # paths to tools
    finemapr_finemap = finemapr_guess_tool("finemap"),
    finemapr_caviar = finemapr_guess_tool("caviar"),
    finemapr_paintor = finemapr_guess_tool("paintor"),
    finemapr_cojo = finemapr_guess_tool("cojo"),
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
