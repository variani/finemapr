
#' Run CAVIAR.
#'
#' @examples
#' ex <- example_finemap()
#' out <- run_caviar(ex$tab1, ex$ld1, args = "-c 2")
#' out <- run_caviar(ex$tab1, ex$ld1, args = "-c 0")
#'
#' @export
run_caviar <- function(tab, ld, 
  dir_run = "run_caviar",
  tool = getOption("finemapr_caviar"), args = "")
{
  ### arg
  stopifnot(class(ld) == "matrix")
  stopifnot(class(dir_run) == "character")
  
  stopifnot(!is.null(rownames(ld)))
  stopifnot(!is.null(colnames(ld)))
  
  stopifnot(!is.null(tool))
  
  ### var
  log_filename <- "finemapr_caviar"
  log_set <- file.path(dir_run, paste0(log_filename, "_set"))
  log_post <- file.path(dir_run, paste0(log_filename, "_post"))
  log_log <- file.path(dir_run, paste0(log_filename, ".log"))
    
  ### process input data: `tab` and `ld`
  tab <- as_data_frame(tab)
  stopifnot(ncol(tab) >= 2)
  names(tab) <- c("snp", "zscore")
 
  tab <- filter(tab, !is.na(zscore)) # exclude missing Z-scores
  
  snps <- tab$snp
  stopifnot(all(snps %in% rownames(ld)))
  stopifnot(all(snps %in% colnames(ld)))
  
  ld <- ld[snps, snps]
  
  stopifnot(all(round(diag(ld), 4) == 1))
  diag(ld) <- 1

  ### create `dir`
  ret_dir_create <- dir.create(dir_run, showWarnings = FALSE, recursive = TRUE)
  
  ### write files
  write_delim(tab, file.path(dir_run, "region.z"), 
    delim = " ", col_names = FALSE)
  write.table(ld, file.path(dir_run, "region.ld"), 
    sep = " ", row.names = FALSE, col.names = FALSE)

  ### run tool
  tool_input <- paste0(args, " -z region.z -l region.ld -o ", log_filename)
  cmd <- paste(tool, tool_input)
  
  dir_cur <- getwd()
  setwd(dir_run)
  
  ret_run <- try({
    system(cmd, input = tool_input)
  })
  
  setwd(dir_cur)
  
  # executed `tool` ok?
  status_run <- ifelse(ret_run == 0, TRUE, FALSE)
  
  # logs
  if(status_run) {
    if(!all(file.exists(log_log, log_post, log_log))) {
      status_run <- FALSE
    }
  }
  
  if(!status_run) {
    out <- list(cmd = cmd, ret = ret_run, status = status_run)
  
    return(out) 
  }
  
  # read output tables
  post <- read_tsv(file.path(dir_run, paste0(log_filename, "_post")))
  
  stopifnot(ncol(post) == 3)
  names(post) <- c("snp", "snp_prob_set", "snp_prob")
  
  post <- arrange(post, -snp_prob)
  
  set <- read_lines(file.path(dir_run, paste0(log_filename, "_set")))

  ### return
  out <- list(cmd = cmd, ret = ret_run, status = status_run, log = log,
    tab = tab, post = post, set = set)
  
  oldClass(out) <- c("FinemaprCaviar", oldClass(out))
   
  return(out) 
}

