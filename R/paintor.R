
#' Run Paintor.
#' 
#' @examples
#' ex <- example_finemap()
#' out <- run_paintor(ex$tab1, ex$ld1)
#'
#' @export
run_paintor <- function(tab, ld, n, annot, annotations, 
  dir_run = "run_paintor",
  tool = getOption("finemapr_paintor"), args = "")
{
  ### arg
  stopifnot(class(ld) == "matrix")
  stopifnot(class(dir_run) == "character")
  
  stopifnot(!is.null(rownames(ld)))
  stopifnot(!is.null(colnames(ld)))
  
  stopifnot(!is.null(tool))
  stopifnot(file.exists(tool))
  
  ### var
  missing_annot <- missing(annot)
  
  missing_n <- missing(n)
  if(!missing_n) {
    num_ind <- n # use another name for `n`
  }
  
  ### process input data: `tab` and `ld`
  tab <- as_data_frame(tab)
  stopifnot(ncol(tab) >= 2)
  names(tab)[c(1, 2)] <- c("snp", "zscore")
 
  tab <- filter(tab, !is.na(zscore)) # exclude missing Z-scores
  
  snps <- tab$snp
  stopifnot(all(snps %in% rownames(ld)))
  stopifnot(all(snps %in% colnames(ld)))
  
  ld <- ld[snps, snps]
  
  stopifnot(all(round(diag(ld), 4) == 1))
  diag(ld) <- 1

  # `annot`
  if(missing_annot) {
    # See comment https://github.com/gkichaev/PAINTOR_V3.0/issues/11#issuecomment-303135031
    # about running PAINTOR without annotations: 
    annot <- data_frame(dummy_ones = rep(1, nrow(tab)))
    
    annotations <- "dummy_ones"
  } else {
    stop("not missing annot")
  }

  ### create `dir`
  ret_dir_create <- dir.create(dir_run, showWarnings = FALSE, recursive = TRUE)

  ### write files
  write_delim(tab, file.path(dir_run, "region"), 
    delim = " ", col_names = TRUE)
  write.table(ld, file.path(dir_run, "region.ld"), 
    sep = " ", row.names = FALSE, col.names = FALSE)
  write_delim(annot, file.path(dir_run, "region.annotations"), 
    delim = " ", col_names = TRUE)
    
  lines_master <- c("region")
  write_lines(lines_master, file.path(dir_run, "region.master"))
  
  ### run tool
  annotations_str <- paste(annotations, collapse = ",")
  tool_input <- paste0("-input region.master -Zhead zscore -LDname ld ",
    "-annotations ", annotations_str, " ", args)
  if(!missing_n) {
    tool_input <- paste0(tool_input, " -num_samples ", num_ind) 
  }
  
  cmd <- paste(tool, tool_input)
  
  dir_cur <- getwd()
  setwd(dir_run)
  
  ret_run <- try({
    system(cmd, input = tool_input)
  })
  
  setwd(dir_cur)
  
  # executed `tool` ok?
  status_run <- ifelse(ret_run == 0, TRUE, FALSE)
  
  # check
  if(!status_run) {
    out <- list(cmd = cmd, ret = ret_run, status = status_run)
  
    return(out) 
  }
  
  #log <- read_lines(log_log)
  
  ### read output tables
  snp <- read_delim(file.path(dir_run, "region.results"), delim = " ")
  
  stopifnot(ncol(snp) == 3)
  names(snp) <- c("snp", "zscore", "snp_prob")
  
  snp <- arrange(snp, -snp_prob) %>%
    mutate(rank = seq(1, n())) %>%
    select(rank, everything())
  
  # enrich
  dat_enrich <- read_delim(file.path(dir_run, "Enrichment.Values"), delim = " ") 
  stopifnot(nrow(dat_enrich) == 1)
  
  gamma_val <- as.numeric(dat_enrich)
  
  # https://github.com/gkichaev/PAINTOR_V3.0/wiki/4.-Interpretation-of-Output#gamma-estimates
  enrich_prob <- c(1 / (1 + exp(gamma_val[1])), 
    1 / (1 + exp(gamma_val[1] + gamma_val[-1])))
  
  enrich <- data_frame(annot = colnames(dat_enrich), 
    gamma = gamma_val, enrich_prob)

  # logBF
  lines_logbf <- read_lines(file.path(dir_run, "Log.BayesFactor"))
  stopifnot(length(lines_logbf) == 1)
  logBF <- as.numeric(lines_logbf)
  
  ### return
  out <- list(cmd = cmd, ret = ret_run, status = status_run, #log = log,
    annotations = annotations, tab = tab, snp = snp, 
    enrich = enrich, logBF = logBF)
  
  oldClass(out) <- c("FinemaprPaintor", "Finemapr", oldClass(out))
   
  return(out) 
}

