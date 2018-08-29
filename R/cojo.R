
#' Run GCTA-COJO.
#' 
#' http://gcta.freeforums.net/thread/178/conditional-joint-analysis-using-summary
#'
#' tab <- "~/git/hemostat/scc/results/sumstats-cojo/sumstats-1.cojo"
#' bed <- "~/git/hemostat/scc/data/FineMapping/forGCTA_onco_r1.bed"
#' 
#' @export
cojo <- function(tab, bed, 
  method = c("select"),
  p = 5e-8,
  dir_run = "run_cojo",
  tool = getOption("finemapr_cojo"), args = "")
{
  ### arg
  method <- match.arg(method)
  stopifnot(class(dir_run) == "character")

  bed <- normalizePath(bed)
  bed <- gsub(".bed$", "", bed)
  
  ### process input data: `tab`
  names_tab <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
  
  tab <- switch(class(tab),
    "character" = read_tsv(tab),
    as_data_frame(tab))
  stopifnot(ncol(tab) == length(names_tab))
  stopifnot(all(names(tab) == names_tab))
  
  snps <- tab$SNP

  ### create `dir`
  ret_dir_create <- dir.create(dir_run, showWarnings = FALSE, recursive = TRUE, mode = "777")
  
  ### write files
  write_tsv(tab, file.path(dir_run, "region.ma"))

  ### run tool
  tool_input <- paste0(args, " --bfile ", bed, " --cojo-file region.ma",
    " --out region", " --cojo-p ", p)
  
  if(method == "select") {
    tool_input <- paste0(tool_input, " --cojo-slct")
  }
  
  cmd <- paste0(tool, tool_input)
  
  dir_cur <- getwd()
  setwd(dir_run)
  
  ret_run <- try({
    system(cmd, input = tool_input)
  })
  
  setwd(dir_cur)
  
  ### read results
  jma <- file.path(dir_run, "region.jma.cojo") %>% read_tsv
  log <- file.path(dir_run, "region.log") %>% read_lines
  badsnps <- file.path(dir_run, "region.freq.badsnps") %>% read_tsv
  
  ### return
  out <- list(cmd = cmd, ret = ret_run, 
    tab = tab,
    jma = jma, log = log, badsnps = badsnps)
  
  oldClass(out) <- c("Cojo", oldClass(out))
  
  return(out) 
}

