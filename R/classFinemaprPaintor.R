#' S3 class FinemaprPaintor.
#'
#' @name FinemaprPaintor
#' @rdname FinemaprPaintor
#'
#' @exportClass FinemaprPaintor

#' @rdname FinemaprPaintor
#' @export
print.FinemaprPaintor <- function(x, ...)
{
  cat(" - cmd:", x$cmd, "\n")
  cat(" - #loci:", x$num_loci, "\n") 
  cat(" - annotations:", paste(x$annotations, collapse = ", "), "\n")
  cat(" - logBF (proportional to the model likelihood):", x$logBF, "\n")
  
  ret <- lapply(seq(1, x$num_loci), function(i) {    
    cat(" - locus:",i, "\n")
    cat("  -- snp:\n")
    print(x$snp[[i]], n = 3)
    cat("  -- ", length(x$snps_credible[[i]]), " snps in ",
      100*x$prop_credible, "% credible set", 
      ": ", paste(x$snps_credible[[i]], collapse = ", "), "...", 
      "\n", sep = "") 
  })
}

#---------------------
# Processing methods
#---------------------

#' @rdname FinemaprPaintor
#' @export
process_annot.FinemaprPaintor <- function(x, annots, annotations, ...)
{
  missing_annot <- missing(annots)
  missing_annotations <- missing(annotations)

  # `annot`
  if(missing_annot) {
    # See comment https://github.com/gkichaev/PAINTOR_V3.0/issues/11#issuecomment-303135031
    # about running PAINTOR without annotations: 
    x$annot <- lapply(seq(1, x$num_loci), function(locus) {
      data_frame(dummy_ones = rep(1, nrow(x$tab[[locus]])))
    })
    x$annotations <- "dummy_ones"
  } else {
    stop("input annotations: not implemented yet")
    
    if(class(annots)[1] != "list") {
      annots <- list(annots)
    }
    stopifnot(length(annots) == x$num_loci)

    out_lds <- lapply(seq_along(annots), function(locus) {
      annot <- annots[[locus]]
    
      # check annotations
      # ...
        
      list(annot = annot)
    })
  

    x$annot <- lapply(out_annot, function(x) x$annot)
    
    stopifnot(!missing_annotations)
    x$annotations <- annotations  
  }

  return(x)
}

#---------------------
# Finemapping methods
#---------------------

#' @rdname FinemaprPaintor
#' @export
write_files.FinemaprPaintor <- function(x, ...)
{
  ### create `dir`
  ret_dir_create <- dir.create(x$dir_run, showWarnings = FALSE, recursive = TRUE)
  #stopifnot(ret_dir_create)

  ### write file of Z-scores
  ret <- lapply(seq_along(x$tab), function(locus) {
    write_delim(
      x$tab[[locus]] %>% select(snp, zscore), 
      file.path(x$dir_run, filename_zscore(x, locus)), 
      delim = " ", col_names = TRUE)
  })
  
  ### write file of ld
  ret <- lapply(seq_along(x$ld), function(locus) {
    write.table(x$ld[[locus]],
      file.path(x$dir_run, filename_ld(x, locus)), 
      sep = " ", row.names = FALSE, col.names = FALSE)
  })

  ### write file of annotations
  ret <- lapply(seq_along(x$annot), function(locus) {
    write_delim(
      x$annot[[locus]],
      file.path(x$dir_run, filename_annot(x, locus)), 
      delim = " ", col_names = TRUE)
  })
  
  ### write master file
  lines_master <- sapply(seq(1, x$num_loci), function(locus) paste0(locus, ".region"))
  write_lines(lines_master, file.path(x$dir_run, filename_master(x)))
}  

#' @rdname FinemaprPaintor
#' @export
run_tool.FinemaprPaintor <- function(x, ...)
{
  annotations_str <- paste(x$annotations, collapse = ",")
  tool_input <- paste0("-input ", filename_master(x),
    " -Zhead ", finemapr_names_tab_zscore(), 
    " -LDname ld ", # hard-coded extension of LD files
    " -annotations ", annotations_str, 
    " ", x$args)
  cmd <- paste(x$tool, tool_input)
  
  dir_cur <- getwd()
  setwd(x$dir_run)
  
  ret_run <- try({
    system(cmd, input = tool_input)
  })
  
  setwd(dir_cur)

  ### return
  x$cmd <- cmd
  x$ret_run <- ret_run
  
  return(x)
}


#' @rdname FinemaprPaintor
#' @export
collect_results.FinemaprPaintor <- function(x, ...)
{
  results <- try({
    lapply(seq(1, x$num_loci), function(locus) {
      #log <- read_lines(file.path(x$dir_run, filename_log(x, locus)))
      
      snp <- file.path(x$dir_run, filename_snp(x, locus)) %>%
          read_delim(, delim = " ", col_types = cols())
      stopifnot(ncol(snp) == 3)
      names(snp) <- c("snp", "zscore", "snp_prob")
      
      snp <- select(snp, snp, snp_prob) %>% 
        arrange(-snp_prob) %>%
        mutate(
          rank_pp = seq(1, n()),
          snp_prob_cumsum = cumsum(snp_prob) / sum(snp_prob)) %>%
        select(rank_pp, snp, snp_prob, snp_prob_cumsum, everything())
      
      snp <- select(x$tab[[locus]], rank_z, snp) %>%
        left_join(snp, ., by = "snp") %>%
        select(rank_z, everything())
      
      list(
        snp = snp)
    })
  })
  
  ### check status and return
  x$status <- ifelse(class(results)[1] == "try-error", 1, 0)
  if(x$status == 0) {
    #x$log <- lapply(results, function(x) x$log)
    x$snp <- lapply(results, function(x) x$snp)
    
    x$snps_credible <- extract_credible_set(x)
  }
  
  return(x)
}

