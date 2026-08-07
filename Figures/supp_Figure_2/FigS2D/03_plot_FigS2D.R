library(gGnome)
library(gTrack)
library(dplyr)

case <- commandArgs(trailingOnly = TRUE)
xG <- case[1] # CCND1, NSD3, EGFR

walks_consAmp_BRCA <- readRDS(paste0('source_data/HLAMPwalks_consensus_ampSum_',xG,'_BRCA_10kb_binned_interpolated.rds'))
walks_consAmp_NSCLC <- readRDS(paste0('source_data/HLAMPwalks_consensus_ampSum_',xG,'_NSCLC_10kb_binned_interpolated.rds'))

CGC <- readRDS('../../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')
gene <- CGC %Q% (Gene_Symbol == xG)
plot_window <- gene + 2e6

#plot
walks_B <- gTrack(walks_consAmp_BRCA,y.field="AmpFraction_smoothed",col="brown3",bars=T)
walks_L <- gTrack(walks_consAmp_NSCLC,y.field="AmpFraction_smoothed",col="royalblue3",bars=T)

# Fisher BRCA vs. NSCLC
sign_BvL_gt <- readRDS(paste0(xG,'_BRCAvsNSCLC_fisher_80consAmp.rds'))

# Fisher left vs. right
sign_B_LvR_gt <- readRDS(paste0(xG,'_BRCA_LvR_fisher_80consAmp.rds'))
sign_L_LvR_gt <- readRDS(paste0(xG,'_NSCLC_LvR_fisher_80consAmp.rds'))

plot(c(gt_ge,
       sign_L_LvR_gt,
       walks_L,
       sign_B_LvR_gt,
       walks_B,
       sign_BvL_gt),plot_window)
