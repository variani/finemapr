#-----------------
# Names in table of Z-scores
#-----------------
finemapr_names_tab_snp <- function() "snp"

finemapr_names_tab_zscore <- function() "zscore"

finemapr_names_tab_pos <- function() "pos"

finemapr_names_tab <- function() 
{
  c(finemapr_names_tab_snp(),
    finemapr_names_tab_zscore())
}

#-----------------
# Find names
#-----------------
finemapr_find_name <- function(target = c("snp", "zscore", "pos"), 
  candidates, strict = FALSE)
{
  target <- match.arg(target)
   
  pat <- switch(target,
    "snp" = "^snp$|^SNP$|^id$|^ID$|^marker|^Marker",
    "zscore" = "^zscore$|^Zscore$|^z$|^Z$",
    "pos" = "^pos|^Pos|^bp$|^BP$",
    stop("error in switch"))
    
  matches <- grep(pat, candidates, value = TRUE)
  if(length(matches) == 0) {
    if(strict) {
      stop("matches for '", target, "' not found; `grep` pattern '", pat, "'")
    }
  }
  if(length(matches) > 1) {
    stop(">1 matches found (`", paste(matches, collapse = ", "), "`); grep pattern '", pat, "'")
  } 
    
  if(length(matches) == 0) {
    matches <- NULL
  }
    
  return(matches)
}
