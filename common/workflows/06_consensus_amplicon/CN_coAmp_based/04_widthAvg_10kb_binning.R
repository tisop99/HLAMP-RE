library(gGnome)
library(data.table)
library(rtracklayer)
library(dplyr)

case <- commandArgs(trailingOnly = TRUE) 
amp <- case[1] # CCND1, NSD3, EGFR
tissue <- case[2] # BRCA, NSCLC

# make a 10kb binning of the reference genome
T2T_seqlengths <- fread('../../../data/refGenome/T2T_contig_lengths_UCSC.tsv',sep='\t')
colnames(T2T_seqlengths) <- c("seqnames","chr_len")
seqlengths <- setNames(T2T_seqlengths$chr_len,T2T_seqlengths$seqnames)
bins <- tileGenome(seqlengths, tilewidth=10000, cut.last.tile.in.chrom=T)

consAmp <- readRDS(paste0('cnBased_consensus_ampSum_extendedCohort_',amp,'_',tissue,'.rds'))

hits <- findOverlaps(bins,consAmp)
ov <- pintersect(bins[queryHits(hits)],consAmp[subjectHits(hits)])
w <- width(ov) * mcols(consAmp)$sample_support[subjectHits(hits)]
signal <- tapply(w, queryHits(hits), sum)

bins$mean_support <- 0
bins$mean_support[as.integer(names(signal))] <- signal / width(bins)[as.integer(names(signal))]
mcols(bins)$AmpFraction <- mcols(bins)$mean_support / max(consAmp$sample_support)

saveRDS(bins,paste0('cnBased_consensus_ampSum_extendedCohort_',amp,'_',tissue,'_10kb_binned.rds'))
