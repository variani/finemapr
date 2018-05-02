
#' Run FINEMAP.
#'
#' @examples
#' ex <- example_finemap()
#' out <- run_finemap(ex$tab1, ex$ld1, ex$n1, args = "--n-causal-max 3")
#' out <- run_finemap(ex$tab1, ex$ld1, ex$n1, args = "--n-causal-max 1")
#'
#' @export
run_finemap <- function(tab, ld, n, 
  dir_run = "run_finemap",
  tool = getOption("finemapr_finemap"), args = "")
{
  ### arg
  stopifnot(class(ld) == "matrix")
  stopifnot(!is.null(rownames(ld)))
  stopifnot(!is.null(colnames(ld)))
  
  ### var
  num_ind <- n # use another name for `n`
  
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

  ### create `dir`
  ret_dir_create <- dir.create(dir_run, showWarnings = FALSE, recursive = TRUE)
  #stopifnot(ret_dir_create)
  
  ### write files
  write_delim(tab[, 1:2], file.path(dir_run, "region.z"), 
    delim = " ", col_names = FALSE)
  write.table(ld, file.path(dir_run, "region.ld"), 
    sep = " ", row.names = FALSE, col.names = FALSE)

  lines_master <- c("z;ld;snp;config;log;n-ind",
    paste0("region.z;region.ld;region.snp;region.config;region.log;", num_ind))
  write_lines(lines_master, file.path(dir_run, "region.master"))
  
  ### run tool
  tool_input <- paste0("--sss --log ", args, " --in-files region.master") 
  cmd <- paste(tool, tool_input)
  
  dir_cur <- getwd()
  setwd(dir_run)
  
  ret_run <- try({
    system(cmd, input = tool_input)
  })
  
  setwd(dir_cur)
  
  # executed `tool` ok?
  status_run <- ifelse(ret_run == 0, TRUE, FALSE)
  
  # read log
  if(status_run) {
    log <- try(read_lines(file.path(dir_run, "region.log")))
    status_run <- ifelse(class(log)[1] == "try-error", FALSE, TRUE)
  }
  
  if(!status_run) {
    out <- list(cmd = cmd, ret = ret_run, status = status_run)
  
    return(out) 
  }
  
  # read output tables
  snp <- read_delim(file.path(dir_run, "region.snp"), delim = " ")
  config <- read_delim(file.path(dir_run, "region.config"), delim = " ")

  # extract output tables
  ncausal <- finemap_extract_ncausal(log)
  
  ### return
  out <- list(cmd = cmd, ret = ret_run, status = status_run, log = log,
    tab = tab, snp = snp, config = config, ncausal = ncausal)
  
  oldClass(out) <- c("FinemaprFinemap", "Finemapr", oldClass(out))
  
  return(out) 
}

finemap_extract_ncausal <- function(log)
{
  lines <- grep("->", log, value = TRUE)
  
  lines <- gsub("\\(|\\)|>", "", lines)
  
  splits <- strsplit(lines, "\\s+")
  
  tab <- data_frame(
    ncausal_num = sapply(splits, function(x) as.integer(x[2])),
    ncausal_prob = sapply(splits, function(x) as.double(x[4])))

  tab <- mutate(tab, type = ifelse(duplicated(ncausal_num), "post", "prior")) 

  return(tab)
}

