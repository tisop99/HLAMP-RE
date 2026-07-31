library(rtracklayer)
library(gGnome)
library(data.table)

case <- commandArgs(trailingOnly = TRUE) # samplesheet file name
amp <- case[1] # CCND1, NSD3, EGFR
tissue <- case[2] # BRCA, NSCLC

if (tissue == "BRCA") {
  tt <- "BRCA"
} else if (tissue == "NSCLC") {
  tt <- "LUSC"
} else {
  stop("Invalid tissue input")
}

### Hartwig/PCAWG
db_gr <- readRDS(paste0('data/joint_HMF_PCAWG_ampSum_cnBased_',amp,'_',tissue,'_T2T.rds'))
max_db <- max(db_gr$sample_support)

# ONT WGS data
np_gr <- readRDS(paste0('ONTsamples_uniStrand_ampSum_CNbased_',amp,'_',tt,'.rds'))
max_np <- max(np_gr$sample_support)

gr_all <- c(db_gr,np_gr)
gr_disjoint <- disjoin(gr_all)
# Find Overlaps between disjoint bins and original ranges (i.e. find origins)
hits <- findOverlaps(gr_disjoint, gr_all)
support_vec <- numeric(length(gr_disjoint))
# Sum the sample_support for each disjoint range 
support_vec <- tapply(
  mcols(gr_all)$sample_support[subjectHits(hits)],
  queryHits(hits),
  sum
)
gr_disjoint$sample_support <- 0
gr_disjoint$sample_support[as.integer(names(support_vec))] <- support_vec

max_both <- max_db + max_np
mcols(gr_disjoint)$AmpFraction <- mcols(gr_disjoint)$sample_support / max_both

saveRDS(gr_disjoint,paste0('cnBased_consensus_ampSum_extendedCohort_',amp,'_',tissue,'.rds'))
