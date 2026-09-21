library(gGnome)
library(gTrack)
library(dplyr)

case <- commandArgs(trailingOnly = TRUE)
xG <- case[1] # CCND1, NSD3, EGFR

CGC <- readRDS('../../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')
gene <- CGC %Q% (Gene_Symbol == xG)
plot_window <- gene + 2e6
gt_ge <- readRDS(paste0('../../../common/data/refGenome/',xG,'_plotting_genetrack_T2T.rds'))

walks_consAmp_BRCA <- readRDS(paste0('source_data/HLAMPwalks_consensus_ampSum_',xG,'_BRCA_10kb_binned.rds'))
walks_consAmp_NSCLC <- readRDS(paste0('source_data/HLAMPwalks_consensus_ampSum_',xG,'_NSCLC_10kb_binned.rds'))

walks_consAmp_BRCA_gr <- dt2gr(walks_consAmp_BRCA)
walks_consAmp_NSCLC_gr <- dt2gr(walks_consAmp_NSCLC)
walks_consAmp_BRCA_gt <- gTrack(walks_consAmp_BRCA_gr, y.field="AmpPercent", col="brown3",bars=T)
walks_consAmp_NSCLC_gt <- gTrack(walks_consAmp_NSCLC_gr, y.field="AmpPercent", col="royalblue3",bars=T)

fisher_bVl <- readRDS(paste0('source_data/',xG,'_BRCAvsNSCLC_FisherStats_CN10_cons80_gr.rds'))
#sign_BvL <- fisher_bVl[fisher_bVl$Fisher_pval<0.05]
sign_BvL <- fisher_bVl[fisher_bVl$qval<0.1]
mcols(sign_BvL)$tissue <- ifelse(mcols(sign_BvL)$BRCA_AmpFreq > mcols(sign_BvL)$NSCLC_AmpFreq, "BRCA", "NSCLC")
colormap_BvL <- list(tissue=c(NSCLC="royalblue3", BRCA="brown3"))
sign_BvL_gt <- gTrack(sign_BvL,gr.colorfield = "color",colormaps = colormap_BvL, height=5)

fisher_LvR_l <- fread(paste0('source_data/',xG,'_NSCLC_LeftVsRight_FisherStats_CN10_cons80.txt'))
fisher_LvR_b <- fread(paste0('source_data/',xG,'_BRCA_LeftVsRight_FisherStats_CN10_cons80.txt'))

fisher_LvR_l_gr <- dt2gr(fisher_LvR_l[, ':='(seqnames = fifelse(right_AmpFreq > left_AmpFreq, right_chr, left_chr),
					                        start = fifelse(right_AmpFreq > left_AmpFreq, right_start, left_start),
					                        end = fifelse(right_AmpFreq > left_AmpFreq, right_end, left_end)
					        		  )])
fisher_LvR_b_gr <- dt2gr(fisher_LvR_b[, ':='(seqnames = fifelse(right_AmpFreq > left_AmpFreq, right_chr, left_chr),
					                        start = fifelse(right_AmpFreq > left_AmpFreq, right_start, left_start),
					                        end = fifelse(right_AmpFreq > left_AmpFreq, right_end, left_end)
								  )])

#sign_LvR_l <- fisher_LvR_l_gr[fisher_LvR_l_gr$Fisher_pval<0.05]
#sign_LvR_b <- fisher_LvR_b_gr[fisher_LvR_b_gr$Fisher_pval<0.05]
sign_LvR_l <- fisher_LvR_l_gr[fisher_LvR_l_gr$qval<0.1]
sign_LvR_b <- fisher_LvR_b_gr[fisher_LvR_b_gr$qval<0.1]
sign_LvR_l_gt <- gTrack(sign_LvR_l, col="royalblue3", height=5)
sign_LvR_b_gt <- gTrack(sign_LvR_b, col="brown3", height=5)

plot(c(gt_ge,
       walks_consAmp_NSCLC_gt,
       sign_LvR_l_gt,
       sign_BvL_gt,
       sign_LvR_b_gt,
       walks_consAmp_BRCA_gt),plot_window)
