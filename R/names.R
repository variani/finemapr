#-----------------
# Names in table of Z-scores
#-----------------
finemapr_names_tab_snp <- function() "snp"

finemapr_names_tab_zscore <- function() "zscore"

finemapr_names_tab <- function() 
{
  c(finemapr_names_tab_snp(),
    finemapr_names_tab_zscore())
}

#-----------------
# Find names
#-----------------
finemapr_find_name <- function(target = c("snp", "zscore"), 
  candidates, strict = FALSE)
{
  target <- match.arg(target)
   
  pat <- switch(target,
    "snp" = "^snp$|^SNP$|^id$|^ID$|^marker$",
    "zscore" = "^zscore$|^Zscore$|^z$|^Z$",
    stop("error in switch"))
    
  matches <- grep(pat, candidates, value = TRUE)
  if(length(matches) == 0) {
    if(strict) {
      stop("matches for '", target, "' not found; `grep` pattern '", pat, "'")
    }
  }
  if(length(matches) > 1) {
    if(strict) {
      stop(">1 matches found (`", paste(matches, collapse = ", "), "`); grep pattern '", pat, "'")
    }
  } 
    
  return(matches)
}
