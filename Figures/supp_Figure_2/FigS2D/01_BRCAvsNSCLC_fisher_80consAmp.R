library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(gGnome)
library(gTrack)

case <- commandArgs(trailingOnly = TRUE)
xG <- case[1] # CCND1

walks_consAmp_BRCA <- readRDS(paste0('source_data/HLAMPwalks_consensus_ampSum_',xG,'_BRCA_10kb_binned_interpolated.rds'))
walks_consAmp_NSCLC <- readRDS(paste0('source_data/HLAMPwalks_consensus_ampSum_',xG,'_NSCLC_10kb_binned_interpolated.rds'))
total_breast <- max(walks_consAmp_BRCA$mean_support)/max(walks_consAmp_BRCA$AmpFraction)
total_lung <- max(walks_consAmp_NSCLC$mean_support)/max(walks_consAmp_NSCLC$AmpFraction)

CGC <- readRDS('../../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')
gene <- CGC %Q% (Gene_Symbol == xG)
target_chr <- as.character(seqnames(gene[1]))

T2T_seqlengths <- fread('../../../common/data/refGenome/T2T_contig_lengths_UCSC.tsv')
colnames(T2T_seqlengths) <- c("seqnames","chr_len")
seqlengths <- setNames(T2T_seqlengths$chr_len,T2T_seqlengths$seqnames)
bins <- tileGenome(seqlengths, tilewidth=10000, cut.last.tile.in.chrom=T)
mcols(bins)$fisher_p <- NA
mcols(bins)$log2FC <- NA
mcols(bins)$fisher_q <- NA

## Fisher test
epsi <- 0.0001 # small pseudocount to avoid division by 0 in foldchange
# define target range (80% consAmp)
b_80 <- walks_consAmp_BRCA[walks_consAmp_BRCA$AmpFraction_smoothed>=0.8]
l_80 <- walks_consAmp_NSCLC[walks_consAmp_NSCLC$AmpFraction_smoothed>=0.8]
start_80 <- min(min(start(b_80)),min(start(l_80)))
end_80 <- max(max(end(b_80)),max(end(l_80)))
target_range <- bins[seqnames(bins)==target_chr]
target_range <- target_range[start(target_range) >= start_80 & end(target_range) <= end_80]
target_range <- keepSeqlevels(target_range, target_chr, pruning.mode = "coarse")
chrs_list <- target_chr

for (chrs in chrs_list) {
  lung <- subsetByOverlaps(walks_consAmp_NSCLC,target_range)
  breast <- subsetByOverlaps(walks_consAmp_BRCA,target_range)
    
  total_lung <- max(lung$mean_support_smoothed)
  total_breast <- max(breast$mean_support_smoothed)
      
  ampCount_lung <- lung$mean_support_smoothed
  ampCount_breast <- breast$mean_support_smoothed
  nonAmpCount_lung <- total_lung - ampCount_lung
  nonAmpCount_breast <- total_breast - ampCount_breast
	  
  fisher_bin <- function(a,b,c,d) {
    mat <- matrix(c(a,b,c,d), nrow = 2)
    ft <- fisher.test(mat) # compute fisher test
    log2fc <- log2(((a+epsi)/total_lung) / ((b+epsi)/total_breast)) # compute foldchange
    c(p.value = ft$p.value, log2FC = log2fc)
  }
	  
  res <- mapply(fisher_bin, ampCount_lung, ampCount_breast, nonAmpCount_lung, nonAmpCount_breast)
  res <- t(res)
  mcols(target_range)$fisher_p <- res[, "p.value"]
  mcols(target_range)$log2FC <- res[, "log2FC"]
  mcols(target_range)$fisher_q <- p.adjust(mcols(target_range)$fisher_p, method = "BH")
}
target_range
min(target_range$fisher_q)
#sign_BvL <- target_range[target_range$fisher_q<0.1]
sign_BvL <- target_range[target_range$fisher_p<0.05]

mcols(sign_BvL)$tissue <- ifelse(mcols(sign_BvL)$log2FC > 0, "NSCLC", "BRCA")
colormap_BvL <- list(tissue=c(NSCLC="royalblue3", BRCA="brown3"))
sign_BvL_gt <- gTrack(sign_BvL,gr.colorfield = "color",colormaps = colormap_BvL, height=5)

saveRDS(sign_BvL_gt,paste0(xG,'_BRCAvsNSCLC_fisher_80consAmp.rds'))
