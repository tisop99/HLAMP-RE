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

CGC <- readRDS('../../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')
gene <- CGC %Q% (Gene_Symbol == xG)
target_chr <- as.character(seqnames(gene[1]))
gene_range <- GRanges(seqnames=target_chr,ranges=IRanges(start=min(start(gene)),end=max(end(gene))))

############################################################################################################################################################ 
## BRCA fisher left-right
b_all <- walks_consAmp_BRCA
total_breast <- max(b_all$mean_support)/max(b_all$AmpFraction)
# get indices of gene within 80% cons Amp
gene_idx <- findOverlaps(gene_range,b_all)
gene_start_idx <- min(subjectHits(gene_idx))
gene_end_idx <- max(subjectHits(gene_idx))

# prepare output
mcols(b_all)$fish_LR_pVal <- NA
mcols(b_all)$fish_LR_log2FC <- 0

# fisher test function
fisher_bin <- function(a,b,c,d) {
  mat <- matrix(c(a,b,c,d), nrow = 2)
  ft <- fisher.test(mat) # compute fisher test
  log2fc <- log2(((a+epsi)/total_breast) / ((b+epsi)/total_breast)) # compute foldchange
  c(p.value = ft$p.value, log2FC = log2fc)
}
# run LvR fisher
for (i in 1:min((gene_start_idx-1),(length(b_all)-gene_end_idx))) {
  if(b_all$AmpFraction_smoothed[gene_start_idx - i] >= 0.8 | b_all$AmpFraction_smoothed[gene_end_idx + i] >= 0.8) {
    ampCount_left <- b_all$mean_support_smoothed[gene_start_idx - i]
    ampCount_right <- b_all$mean_support_smoothed[gene_end_idx + i]
    nonAmpCount_left <- total_breast - ampCount_left
    nonAmpCount_right <- total_breast - ampCount_right
	    
    res <- fisher_bin(ampCount_left, ampCount_right, nonAmpCount_left, nonAmpCount_right)
    res <- t(res)
    if(res[, "log2FC"] > 0 ) {
      mcols(b_all)$fish_LR_pVal[gene_start_idx - i] <- res[, "p.value"]
      mcols(b_all)$fish_LR_log2FC[gene_start_idx - i] <- res[, "log2FC"]
    }
    if(res[, "log2FC"] < 0 ) {
      mcols(b_all)$fish_LR_pVal[gene_end_idx + i] <- res[, "p.value"]
      mcols(b_all)$fish_LR_log2FC[gene_end_idx + i] <- res[, "log2FC"]
    } 
  }
}
# compute q-value
b_LvR <- b_all[!is.na(b_all$fish_LR_pVal)]
mcols(b_LvR)$fish_LR_qVal <- 1
mcols(b_LvR)$fish_LR_qVal <- p.adjust(b_LvR$fish_LR_pVal, method = "BH")

saveRDS(b_LvR,paste0(xG,'_BRCA_LvR_fisher_80consAmp.rds'))
############################################################################################################################################################ 
## NSCLC fisher left-right
l_all <- walks_consAmp_NSCLC
total_lung <- max(l_all$mean_support)/max(l_all$AmpFraction)
# get indices of gene within 80% cons Amp
gene_idx <- findOverlaps(gene_range,l_all)
gene_start_idx <- min(subjectHits(gene_idx))
gene_end_idx <- max(subjectHits(gene_idx))

# prepare output
mcols(l_all)$fish_LR_pVal <- NA
mcols(l_all)$fish_LR_log2FC <- 0

# fisher test function
fisher_bin <- function(a,b,c,d) {
	  mat <- matrix(c(a,b,c,d), nrow = 2)
  ft <- fisher.test(mat) # compute fisher test
    log2fc <- log2(((a+epsi)/total_lung) / ((b+epsi)/total_lung)) # compute foldchange
    c(p.value = ft$p.value, log2FC = log2fc)
}

# run LvR fisher
for (i in 1:min((gene_start_idx-1),(length(l_all)-gene_end_idx))) {
  if(l_all$AmpFraction_smoothed[gene_start_idx - i] >= 0.8 | l_all$AmpFraction_smoothed[gene_end_idx + i] >= 0.8) {
    ampCount_left <- l_all$mean_support_smoothed[gene_start_idx - i]
    ampCount_right <- l_all$mean_support_smoothed[gene_end_idx + i]
    nonAmpCount_left <- total_lung - ampCount_left
    nonAmpCount_right <- total_lung - ampCount_right
	    
    res <- fisher_bin(ampCount_left, ampCount_right, nonAmpCount_left, nonAmpCount_right)
    res <- t(res)
    if(res[, "log2FC"] > 0 ) {
      mcols(l_all)$fish_LR_pVal[gene_start_idx - i] <- res[, "p.value"]
      mcols(l_all)$fish_LR_log2FC[gene_start_idx - i] <- res[, "log2FC"]
    }
    if(res[, "log2FC"] < 0 ) {
      mcols(l_all)$fish_LR_pVal[gene_end_idx + i] <- res[, "p.value"]
      mcols(l_all)$fish_LR_log2FC[gene_end_idx + i] <- res[, "log2FC"]
    } 
  }
}
# compute q-value
l_LvR <- l_all[!is.na(l_all$fish_LR_pVal)]
mcols(l_LvR)$fish_LR_qVal <- 1
mcols(l_LvR)$fish_LR_qVal <- p.adjust(l_LvR$fish_LR_pVal, method = "BH")

saveRDS(l_LvR,paste0(xG,'_NSCLC_LvR_fisher_80consAmp.rds'))
############################################################################################################################################################ 
