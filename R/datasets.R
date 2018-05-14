#---------------------------
# Datasets from finemap
#---------------------------

#' @export
example_finemap <- function(dir_example = "~/apps/finemap/example/")
{
  ### read
  master <- read_delim(file.path(dir_example, "data"), 
    delim = ";", col_names = TRUE, col_types = cols())

  ### dataset 1
  tab1 <- read_delim(file.path(dir_example, "region1.z"), 
    delim = " ", col_names = FALSE, col_types = cols()) 
  names(tab1) <- c("snp", "zscore")  
  
  ld1 <- read_delim(file.path(dir_example, "region1.ld"), 
    delim = " ", col_names = FALSE, col_types = cols()) 
  ld1 <- as.matrix(ld1)
  rownames(ld1) <- tab1$snp
  colnames(ld1) <- tab1$snp
  
  n1 <- master[["n-ind"]][1]

  ### dataset 2
  tab2 <- read_delim(file.path(dir_example, "region2.z"), 
    delim = " ", col_names = FALSE, col_types = cols()) 
  names(tab2) <- c("snp", "zscore")  
  
  ld2 <- read_delim(file.path(dir_example, "region2.ld"), 
    delim = " ", col_names = FALSE, col_types = cols()) 
  ld2 <- as.matrix(ld2)
  rownames(ld2) <- tab2$snp
  colnames(ld2) <- tab2$snp
  
  n2 <- master[["n-ind"]][2]
    
  ### return
  list(tab1 = tab1, ld1 = ld1, n1 = n1,
    tab2 = tab2, ld2 = ld2, n2 = n2)
}
