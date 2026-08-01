library(gGnome)
library(data.table)
library(rtracklayer)
library(dplyr)

case <- commandArgs(trailingOnly = TRUE) 
amp <- case[1] # CCND1, NSD3, EGFR
tissue <- case[2] # BRCA / NSCLC

# create a 10kb binning from ref genome
T2T_seqlengths <- fread('/charite-store-f/f-cc05-nanoamps/T2T/T2T_contig_lengths_USCS.tsv',sep='\t')
colnames(T2T_seqlengths) <- c("seqnames","chr_len")
seqlengths <- setNames(T2T_seqlengths$chr_len,T2T_seqlengths$seqnames)
bins <- tileGenome(seqlengths, tilewidth=10000, cut.last.tile.in.chrom=T)

# load the consensus amplicon
consAmp <- readRDS(paste0('HLAMPwalks_consensus_ampSum_',amp,'_',tissue,'.rds'))

hits <- findOverlaps(bins,consAmp)
ov <- pintersect(bins[queryHits(hits)],consAmp[subjectHits(hits)])
w <- width(ov) * mcols(consAmp)$sample_support[subjectHits(hits)]
signal <- tapply(w, queryHits(hits), sum)

# compute mean sample_support and amp fraction
bins$mean_support <- 0
bins$mean_support[as.integer(names(signal))] <- signal / width(bins)[as.integer(names(signal))]
mcols(bins)$AmpFraction <- mcols(bins)$mean_support / max(consAmp$sample_support)

saveRDS(bins,paste0('HLAMPwalks_consensus_ampSum_',amp,'_',tissue,'_10kb_binned.rds'))
