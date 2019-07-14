#' S3 class FinemaprFinemap.
#'
#' @name FinemaprFinemap
#' @rdname FinemaprFinemap
#'
#' @exportClass FinemaprFinemap

#---------------------
# Finemapping methods
#---------------------

#' @rdname FinemaprFinemap
#' @export
write_files.FinemaprFinemap <- function(x, ...)
{
  ### create `dir`
  ret_dir_create <- dir.create(x$dir_run, showWarnings = FALSE, recursive = TRUE)
  #stopifnot(ret_dir_create)

  ### write file of Z-scores
  ret <- lapply(seq_along(x$tab)[[1]], function(locus) {
    write_delim(
      prepare_zscore_writing(x$tab[[2]]), # zscore
      file.path(x$dir_run, filename_zscore(x, locus)), 
      delim = " ", col_names = FALSE)
  })
  
  ### write file of ld
  ret <- lapply(seq_along(x$ld), function(locus) {
    write.table(x$ld[[locus]],
      file.path(x$dir_run, filename_ld(x, locus)), 
      sep = " ", row.names = FALSE, col.names = FALSE)
  })

  ### write master file
  lines_master <- c(
    paste0("z;ld;snp;config",
      ifelse(is.null(x$prior_k), "", ";k"),
      ";log;n-ind"),
    sapply(seq(1, x$num_loci), function(locus) {
      paste0(
        filename_zscore(x, locus), ";", 
        filename_ld(x, locus), ";",
        filename_snp(x, locus), ";",
        filename_config(x, locus), ";",
        ifelse(is.null(x$prior_k), "", paste0(filename_k(x, locus), ";")),
        filename_log(x, locus), ";",
        x$n[[locus]])
    }))
  write_lines(lines_master, file.path(x$dir_run, filename_master(x)))
  
  ### write optional files
  if(!is.null(x$prior_k)) {
    ret <- lapply(seq(1, x$num_loci), function(locus) {
      write_lines(paste(x$prior_k, collapse = " "), file.path(x$dir_run, filename_k(x, locus)))
    })
  }
}  

#' @rdname FinemaprFinemap
#' @export
run_tool.FinemaprFinemap <- function(x, ...)
{
  tool_input <- paste0("--sss --log ", 
    ifelse(is.null(x$prior_k), "", " --prior-k "),
    x$args, " --in-files ", filename_master(x))
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

#' @rdname FinemaprFinemap
#' @export
collect_results.FinemaprFinemap <- function(x, ...)
{
  results <- try({
    lapply(seq(1, x$num_loci), function(locus) {
      log <- read_lines(file.path(x$dir_run, sub("[0-9]*\\.", "", filename_log(x, locus))))
      
      snp <- file.path(x$dir_run, sub("[0-9]*\\.", "", filename_snp(x, locus))) %>%
          read_delim(, delim = " ", col_types = cols())
      snp <- snp[locus,]
      
      snp <- arrange(data.table(snp), as.numeric(-snp_prob)) %>%
        mutate(
          rank_pp = seq(1, n()),
          snp_prob_cumsum = cumsum(snp_prob) / sum(snp_prob)) %>%
        select(rank_pp, snp, snp_prob, snp_prob_cumsum, snp_log10bf) #, snp_log10bf)
      

      snp <- merge(data.table(x$tab[which(x$tab[[1]]==snp$snp),]), snp)

      config_list <- file.path(x$dir_run, sub("[0-9]*\\.", "", filename_config(x, locus))) %>%
          read_delim(, delim = " ", col_types = cols())
      
      list(
        log = log,
        snp = snp,
        config = config_list[locus,],
        ncausal = finemap_extract_ncausal(log))
    })
  })
  
  ### check status and return
  x$status <- ifelse(class(results)[1] == "try-error", 1, 0)
  if(x$status == 0) {
    x$log <- lapply(results, function(x) x$log)
    x$snp <- lapply(results, function(x) x$snp)
    x$config <- lapply(results, function(x) x$config)
    x$ncausal <- lapply(results, function(x) x$ncausal)
    
    x$snps_credible <- extract_credible_set(x)
  }
  
  return(x)
}

#---------------------
# Print/plot methods
#---------------------

#' @rdname FinemaprFinemap
#' @export
print.FinemaprFinemap <- function(x, ...)
{
  cat(" - tables of results: `config`, `snp`, `ncausal`\n")
  
  # ret <- lapply(seq(1, x$num_loci), function(i) {    
  #   cat(" - locus:",i, "\n")
  #   cat("  -- config:\n")
  #   cat("  -- input snps: ", length(x$snps_finemap[[i]]), " fine-mapped",
  #     " + ", length(x$snps_missing_finemap[[i]]), " missing Z/LD",
  #     " = ", length(x$snps_zscore[[i]]), " in total\n", sep = "")
  #   print(x$config, n = 3)
  #   cat("  -- snp:\n")
  #   print(x$snp[[2]][i])
  #   cat("  -- ", length(x$snps_credible[[i]]), " snps in ",
  #     100*x$prop_credible, "% credible set", 
  #     ": ", paste(x$snps_credible[[i]], collapse = ", "), "...", 
  #     "\n", sep = "") 
  # })

  # return(invisible())
  
  cat(" - command:", x$cmd, "\n")
    
  if(x$status) {
    cat(" - see log output in `log`\n")
    cat(" - tables of results: `config`, `snp`, `ncausal`\n")
    
    cat(" - config:\n")
    print(x$config, n = 3)

  }
}

#' @rdname FinemaprFinemap
#' @export
plot.FinemaprFinemap <- function(x, 
  grid_nrow = NULL, grid_ncol = NULL, 
  label_size = getOption("finemapr_label_size"),
  label_size_config = label_size, label_size_snp = label_size,  
  top_rank = getOption("top_rank"),
  top_rank_config = top_rank, top_rank_snp = top_rank,   
  lim_prob = getOption("lim_prob"),
  lim_prob_config = lim_prob, lim_prob_snp = lim_prob, lim_prob_ncausal = lim_prob,
  ...)
{
  p1 <- plot_ncausal(x, 
    lim_prob = lim_prob_ncausal, ...)
  p2 <- plot_config(x,  
    top_rank = top_rank_config, 
    label_size = label_size_config, 
    lim_prob = lim_prob_config, ...)
  p3 <- plot_snp(x, 
    top_rank = top_rank_snp,
    label_size = label_size_snp, 
    lim_prob = lim_prob_snp, ...)
  
  plot_grid(p1, p2, p3, nrow = grid_nrow, ncol = grid_ncol)
}

#---------------
# Custom methods
#---------------

#' @rdname FinemaprFinemap
#' @export
plot_ncausal.FinemaprFinemap <- function(x, locus = 1,
  lim_prob = c(0, 1), # automatic limits
  ...)
{
  ptab <- x$ncausal
  
  sum_prop_zero <- filter(ptab, ncausal_num == 0)[["prob"]]  %>% sum
  if(sum_prop_zero == 0) {
    ptab <- filter(ptab, ncausal_num != 0)
  }
  
  ptab <- mutate(ptab, 
    ncausal_num = factor(ncausal_num, levels = sort(unique(ncausal_num), decreasing = TRUE)),
    type = factor(type, levels = c("prior", "post")))
    
  p <- ggplot(ptab, aes(ncausal_num, ncausal_prob, fill = type)) + 
    geom_hline(yintercept = 1, linetype = 3) + 
    geom_bar(stat = "identity", position = "dodge") + 
    coord_flip() + theme(legend.position = "top") + 
    scale_fill_manual(values = c("grey50", "orange")) +
    ylim(lim_prob)
    
  return(p)
}

#' @rdname FinemaprFinemap
#' @export
plot_config.FinemaprFinemap <- function(x, locus = 1,
  lim_prob = c(0, 1.5), 
  label_size = getOption("finemapr_label_size"),  
  top_rank = getOption("top_rank"),  
  ...)
{
  ptab <- x$config

  ptab <- head(ptab, top_rank)

  ptab <- mutate(ptab,
    label = paste0(config, "\n", 
      "P = ", round(config_prob, 2),
      "; ", "log10(BF) = ", round(config_log10bf, 2)))

  ggplot(ptab, aes(config_prob, rank)) + 
    geom_vline(xintercept = 1, linetype = 3) + 
    geom_point() + 
    geom_segment(aes(xend = config_prob, yend = rank, x = 0)) + 
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.025, size = label_size) + 
    xlim(lim_prob) + 
    scale_y_continuous(limits  = c(top_rank + 0.5, 0.5), trans = "reverse")
}

#' @rdname FinemaprFinemap
#' @export
plot_snp.FinemaprFinemap <- function(x, locus = 1,
  lim_prob = c(0, 1.5), 
  label_size = getOption("finemapr_label_size"),  
  top_rank = getOption("top_rank"),  
  ...)
{
  ptab <- x$snp

  ptab <- head(ptab, top_rank)

  ptab <- mutate(ptab,
    rank = seq(1, n()), 
    label = paste0(snp, "\n", 
      "P = ", round(snp_prob, 2),
      "; ", "log10(BF) = ", round(snp_log10bf, 2)))

  ggplot(ptab, aes(snp_prob, rank)) +
    geom_vline(xintercept = 1, linetype = 3) + 
    geom_point() + 
    geom_segment(aes(xend = snp_prob, yend = rank, x = 0)) + 
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.025, size = label_size) + 
    xlim(lim_prob) + 
    scale_y_continuous(limits  = c(top_rank + 0.5, 0.5), trans = "reverse")
}
