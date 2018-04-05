#'
#'
#' https://data.broadinstitute.org/alkesgroup/sumstats_formatted/readme.txt
#' https://github.com/bulik/ldsc/wiki/Partitioned-Heritability
#' 
#' @export
sumstats_bmi <- function()
{
  url <- "https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_BMI1.sumstats"

  read_tsv(url, col_types = cols()) %>%
    mutate(P = pchisq(CHISQ, df = 1, lower.tail = FALSE)) %>%
    select(-CHISQ)
}
