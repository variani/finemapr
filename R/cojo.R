
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
  
  tab <- switch(class(tab)[1],
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
  
  snps_index <- jma$SNP
  
  ### run conditional analysis for each index snps (snps_index)
  cond <- lapply(seq_along(snps_index), function(i) {
    snp_i <- snps_index[i]
    snps_cond <- snps_index[-i]
    out_i <- paste0("region_", i)
    snplist_i <- paste0("cond.snplist_", i)
    
    write_lines(snps_cond, file.path(dir_run, snplist_i ))
    
    tool_input <- paste0(args, " --bfile ", bed, " --cojo-file region.ma",
      paste0(" --out ", out_i), paste0(" --cojo-cond ", snplist_i)) #--cojo-collinear 0.99")
  
    cmd <- paste0(tool, tool_input)
  
    dir_cur <- getwd()
    setwd(dir_run)
  
    ret_run <- try({
      system(cmd, input = tool_input)
    })
    
    setwd(dir_cur)
    
    # read results
    list(
      snp_index = snp_i, snps_cond = snps_cond,
      cma = read_tsv(file.path(dir_run, paste0(out_i, ".cma.cojo"))))
  })

  ### compute credible set using ABF/FINEMAP 1
  cond <- lapply(cond, function(x) {
    tab <- x$cma
    tab <- mutate(tab, zscore = bC/bC_se)
   
    ld <- simulate_ld_diag(tab$SNP) 
    fm <- finemapr(tab, ld, round(mean(tab$n), 0), args = "--n-causal-max 1")
    
    #c(x, list(bf = fm$snp[[1]], snps_credible = fm$snps_credible[[1]]))
    c(x, list(fm = fm))
  })
    
  ### return
  out <- list(cmd = cmd, ret = ret_run, 
    tab = tab,
    jma = jma, log = log, badsnps = badsnps,
    snps = snps, snps_index = snps_index,
    cond = cond)
  
  oldClass(out) <- c("Cojo", oldClass(out))
  
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

