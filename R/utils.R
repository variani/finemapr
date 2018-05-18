
prepare_zscore_writing <- function(tab)
{
  filter(tab, finemap) %>% select(snp, zscore)
}

merge_tab_snp <- function(tab, snp)
{
    select(tab, rank_z, snp) %>%
    left_join(snp, by = "snp") %>%
    select(snp, rank_z, rank_pp, everything()) %>%
    arrange(rank_pp)
}
