##
# genome data from PCAWG,HMF were hg19 based and lifted to T2T
##
library(gGnome)
library(data.table)
library(rtracklayer)
library(zoo)

case <- commandArgs(trailingOnly = TRUE)
amp <- case[1] # CCND1, NSD3, EGFR
tissue <- case[2] # BRCA, NSCLC

if (tissue == "NSCLC") {
  tt <- "LUSC"
} else if (tissue == "BRCA") {
  tt <- "BRCA"
} else {
  print("invalid tissue type supplied")
}

consAmp <- readRDS(paste0('cnBased_consensus_ampSum_extendedCohort_',amp,'_',tissue,'_10kb_binned.rds'))
db <- readRDS(paste0('data/joint_HMF_PCAWG_ampSum_cnBased_',amp,'_',tissue,'_T2T.rds'))
ont <- readRDS(paste0('ONTsamples_uniStrand_ampSum_CNbased_',amp,'_',tt,'.rds'))

### total liftOver holes: 0 coverage from lifted db, but coverage from ont samples (fully covered)
lifthole <- which(countOverlaps(consAmp,db) == 0 & (consAmp %within% ont))
mcols(consAmp)$lifthole <- FALSE
mcols(consAmp)$lifthole[lifthole] <- TRUE
lifthole_logical_vec <- mcols(consAmp)$lifthole  ## logical vector for later
###
### partial liftOver holes: 100% covered by ONT samples, but only partially covered by lifted db samples ++
### for which there is a local depression, i.e. observed mean_support drops by atleast [[ 30% ]] compared to the expected support (roll_median in 250kb window) ++
### for which the expected signal is not consistent with the coverage (i.e. cannot be explained by the coverage drop of the lifted samples) 
observed <- mcols(consAmp)$mean_support

k <- 25 ## 250kb rolling window (10kb bins)
local_signal <- rollmedian(observed, k=k, fill=NA, align="center")
mcols(consAmp)$local_signal <- local_signal

# get bin coverage of PCAWG/HMF db samples
hits_lifted <- findOverlaps(consAmp, db, ignore.strand=T)
lifted_bp <- width(pintersect(consAmp[queryHits(hits_lifted)],db[subjectHits(hits_lifted)]))
lifted_bp_sum <- tapply(lifted_bp, queryHits(hits_lifted), sum)
lifted_bp_vec <- numeric(length(consAmp))
lifted_bp_vec[as.integer(names(lifted_bp_sum))] <- lifted_bp_sum
mcols(consAmp)$lifted_bp_frac <- lifted_bp_vec / 10000
coverage_lifted <- mcols(consAmp)$lifted_bp_frac
# get bin coverage of ONT samples
hits_ont <- findOverlaps(consAmp, ont, ignore.strand=T)
ont_bp <- width(pintersect(consAmp[queryHits(hits_ont)],ont[subjectHits(hits_ont)]))
ont_bp_sum <- tapply(ont_bp, queryHits(hits_ont), sum)
ont_bp_vec <- numeric(length(consAmp))
ont_bp_vec[as.integer(names(ont_bp_sum))] <- ont_bp_sum
mcols(consAmp)$ont_bp_frac <- ont_bp_vec / 10000
coverage_ont <- mcols(consAmp)$ont_bp_frac
# condition: bin is fully covered by ONT samples & less than 90% covered by lifted db samples (& exclude 0 case since already covered above)
coverage_incomplete <- coverage_ont==1 & coverage_lifted < 0.9 & coverage_lifted > 0

# define dips between observed and expected support
drop_ratio <- observed / (local_signal + 1e-6)  # add small count to avoid division by 0
dip <- drop_ratio * 0.7 ## drop by atleast 30%

# expected support given the coverage
expected <- local_signal * coverage_lifted

exp <- observed / (expected + 1e-6)  # add small count to avoid division by 0
exp_by_cov <- abs(exp - 1) < 0.15 ## 0.15 is heuristik
mcols(consAmp)$exp_by_cov <- exp_by_cov

# get indices for which these conditions are true
idx <- which(lifthole_logical_vec | (coverage_incomplete & dip & exp_by_cov))
mcols(consAmp)$hole <- F
mcols(consAmp)$hole[idx] <- T

###
### interpolate identified liftOver holes with mean of second nearest non-hole neighbors (the nearest non-hole neighbors might also be affected)
mcols(consAmp)$mean_support_smoothed <- mcols(consAmp)$mean_support
gr1_list <- split(consAmp, seqnames(consAmp)) # do this for each chromosome separately

gr1_list <- lapply(gr1_list, function(gr) {
  no_idx <- which(gr$hole == FALSE)
  yes_idx <- which(gr$hole == TRUE)

  pos_no <- start(gr)[no_idx]
  val_no <- gr$mean_support_smoothed[no_idx]
  pos_yes <- start(gr)[yes_idx]

  # find closes non hole:
  i <- findInterval(pos_yes,pos_no)
  # second nearest neighbor since nearest neighbor also affected
  n <- length(pos_no)
  left2 <- i-1
  right2 <- i+2
  valid <- left2 >= 1 & right2 <=n

  new_vals <- rep(NA_real_, length(yes_idx))
  new_vals[valid] <- (val_no[left2[valid]] + val_no[right2[valid]]) / 2
  gr$mean_support_smoothed[yes_idx] <- new_vals

  # also overwrite nearest neighbor
  gr$mean_support_smoothed[no_idx[i[valid]]] <- new_vals[valid]
  gr$mean_support_smoothed[no_idx[i[valid] +1]] <- new_vals[valid]

  gr
})

consAmp <- unlist(GRangesList(gr1_list), use.names=F)
max_N <- (max(mcols(consAmp)$mean_support) / max(mcols(consAmp)$AmpFraction))
mcols(consAmp)$AmpFraction_smoothed <- mcols(consAmp)$mean_support_smoothed / max_N

saveRDS(consAmp,paste0('cnBased_consensus_ampSum_extendedCohort_',amp,'_',tissue,'_10kb_binned_interpolated.rds'))
