
#' Run GCTA-COJO.
#' 
#' http://gcta.freeforums.net/thread/178/conditional-joint-analysis-using-summary
#'
#' tab <- "~/git/hemostat/scc/results/sumstats-cojo/sumstats-1.cojo"
#' bed <- "~/git/hemostat/scc/data/FineMapping/forGCTA_onco_r1.bed"
#' 
#' @export
cojo <- function(tab, bed, 
  method = c("select", "cond"),
  snps_cond = NULL,
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
  
  tab <- switch(class(tab)[1],
    "character" = read_tsv(tab, col_types = "cccnnnnn"),
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
    " --out region")
  
  if(method == "select") {
    tool_input <- paste0(tool_input, " --cojo-slct")
  } else if(method == "cond") {
    stopifnot(!is.null(snps_cond))
    write_lines(snps_cond, file.path(dir_run, "cond.snplist"))
    
    tool_input <- paste0(tool_input, " --cojo-cond cond.snplist")
  }
  
  cmd <- paste0(tool, tool_input)
  
  dir_cur <- getwd()
  setwd(dir_run)
  
  ret_run <- try({
    system(cmd, input = tool_input)
  })
  
  setwd(dir_cur)
  
  ### read results
  log <- file.path(dir_run, "region.log") %>% read_lines
  badsnps <- file.path(dir_run, "region.badsnps") %>% read_tsv(col_types = "cccc")
  badfreqs <- file.path(dir_run, "region.freq.badsnps") %>% read_tsv(col_types = "ccccnn")
    
  jma <- snps_index <- cma <- NULL
  if(method == "select") {
    jma <- file.path(dir_run, "region.jma.cojo") %>% read_tsv(col_types = "ccncnnnnnnnnnn")
    snps_index <- jma$SNP
  } else if(method == "cond") {
    cma <- read_tsv(file.path(dir_run, "region.cma.cojo"))
  }   
  
  ### return
  out <- list(cmd = cmd, ret = ret_run, 
    tab = tab,
    # select
    jma = jma, log = log, badsnps = badsnps, badfreqs = badfreqs,
    snps = snps, snps_index = snps_index,
    # cond
    snps_cond = snps_cond, cma = cmd)
  
  oldClass(out) <- c("Cojo", oldClass(out))
  
  return(out) 
}

#' @export
run_cojo <- function(tab, bed,
  args = "", args2 = "", ...)
{
  # step 1: select index snps
  cojo_select <- cojo(tab, bed,  method = "select", args = args, ...)
  
  snps_index <- cojo_select$snps_index
  
  ### step 2: conditional analysis for each index snps (`snps_index`)
  cond <- lapply(seq_along(snps_index), function(i) {
    snp_i <- snps_index[i]
    snps_cond <- snps_index[-i]
    
    if(length(snps_cond)) {
      stop("not implemented")
    } else {
      cma <- cojo_select$tab
    }
    
    # abf
    abf <- abf(cma$b, cma$se, cma$SNP)
    
    snp_below <- abf %>% filter(snp_prob_cumsum <= 0.95)
    snps_credible <- head(abf, nrow(snp_below) + 1) %$% snp
    
    # read results
    list(
      snp_index = snp_i, snps_cond = snps_cond,
      cma = cma, abf = abf, snps_credible = snps_credible)
  })

  ### return
  out <- cojo_select
  out$cond <- cond
 
  return(out) 
}

plot.Cojo <- function(x, locus = 1, digits = 1)
{
  snp_index <- x$cond[[locus]]$snp_index
  p <- subset(x$jma, SNP == snp_index , select = "p", drop = TRUE)
  pJ <- subset(x$jma, SNP == snp_index , select = "pJ", drop = TRUE)
  pC <- subset(x$cond[[locus]]$cma, SNP == snp_index , select = "pC", drop = TRUE)
  
  str_index_credible <- ifelse(
    snp_index %in% x$cond[[locus]]$fm$snps_credible,
    "(inside credible set)",
    "(outside credible set)")
  
  str_pval <- paste0("p = ", format.pval(p, digits = digits),
    "; pJ = ", format.pval(pJ, digits = digits), "; ", 
    "pC = ", format.pval(pC, digits = digits))
  
  title <- paste0("Index SNP #", locus, ": ", x$cond[[locus]]$snp_index,
    " ", str_index_credible)
  
  plot_zscore(x$cond[[locus]]$fm, selected = snp_index) + 
    labs(title = title, subtitle = str_pval)
}

