
#' @export
run_finemap <- function(ztab, ld, dir_run = "finemap")
{
  ### arg
  stopifnot(class(ld) == "matrix")
  stopifnot(!is.null(rownames(ld)))
  stopifnot(!is.null(colnames(ld)))
  
  ### var
  ztab <- as_data_frame(ztab)
  stopifnot(ncol(ztab) == 2)
  names(ztab) <- c("snp", "zscore")
  
  ### process input data: `ztab` and `ld`
  ztab <- filter(ztab, !is.na(ztab)) # exclude missing Z-scores
  
  snps <- ztab$snp
  stopifnot(all(snps %in% rownames(ld)))
  stopifnot(all(snps %in% colnames(ld)))
  
  ld <- ld[snps, snps]
  
  stopifnot(all(round(diag(ld), 4) == 1))
  diag(ld) <- 1

  ### create `dir`
  ret_dir_create <- dir.create(dir_run, showWarnings = FALSE)
  stopifnot(ret_dir_create)
  
  ### write files
  write_delim(ztab, "region.z", delim = " ", col_names = FALSE)
  write.table(ld, "region.ld", sep = " ", row.names = FALSE, col.names = FALSE)

  lines_master <- c("z;ld;snp;config;log;n-ind",
    "region.z;region.ld;region.snp;region.config;region.log;3300")
  write_lines(lines_master, "region.master")
  
  ### run tool
  #ret_run_tool <- 
}
