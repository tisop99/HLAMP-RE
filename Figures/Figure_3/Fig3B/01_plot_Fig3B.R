library(GxG)
library(gGnome)
library(rtracklayer)
library(data.table)

case <- commandArgs(trailingOnly = TRUE)
xG <- case[1] # CCND1, NSD3, EGFR

CGC <- readRDS('../../../common/data/refGenome/cgc_T2T_gr_UCSC.rds')
gene <- CGC %Q% (Gene_Symbol == xG)
plot_window <- gene + 10e5
gt_ge <- readRDS(paste0('../../../common/data/refGenome/',xG,'_plotting_genetrack_T2T.rds'))

# consensus amp BRCA/NSCLC
walks_consAmp_BRCA <- readRDS(paste0('source_data/HLAMPwalks_consensus_ampSum_',xG,'_BRCA_10kb_binned.rds'))
walks_consAmp_NSCLC <- readRDS(paste0('source_data/HLAMPwalks_consensus_ampSum_',xG,'_NSCLC_10kb_binned.rds'))
walks_consAmp_BRCA_gr <- dt2gr(walks_consAmp_BRCA)
walks_consAmp_NSCLC_gr <- dt2gr(walks_consAmp_NSCLC)
walks_consAmp_BRCA_gt <- gTrack(walks_consAmp_BRCA_gr, y.field="AmpPercent", col="brown3")
walks_consAmp_NSCLC_gt <- gTrack(walks_consAmp_NSCLC_gr, y.field="AmpPercent", col="royalblue3")

# HLAMP asymmetries (FisherStats)
fisher_bVl <- readRDS(paste0('source_data/',xG,'_BRCAvsNSCLC_FisherStats_CN10_cons80_gr.rds'))
sign_BvL <- fisher_bVl[fisher_bVl$qval<0.1]
mcols(sign_BvL)$tissue <- ifelse(mcols(sign_BvL)$BRCA_AmpFreq > mcols(sign_BvL)$NSCLC_AmpFreq, "BRCA", "NSCLC")
colormap_BvL <- list(tissue=c(NSCLC="royalblue3", BRCA="brown3"))
sign_BvL_gt <- gTrack(sign_BvL,gr.colorfield = "color",colormaps = colormap_BvL, height=5)

fisher_LvR_b <- fread(paste0('source_data/',xG,'_BRCA_LeftVsRight_FisherStats_CN10_cons80.txt'))
fisher_LvR_b_gr <- dt2gr(fisher_LvR_b[, ':='(seqnames = fifelse(right_AmpFreq > left_AmpFreq, right_chr, left_chr),
                                                                start = fifelse(right_AmpFreq > left_AmpFreq, right_start, left_start),
                                                                end = fifelse(right_AmpFreq > left_AmpFreq, right_end, left_end)
                                                                )])
sign_LvR_b <- fisher_LvR_b_gr[fisher_LvR_b_gr$qval<0.1]
sign_LvR_b_gt <- gTrack(sign_LvR_b, col="brown3", height=5)

# DSS differential analysis
dmrs4 <- readRDS(paste0('source_data/',xG,'_DSS_BRCAvsLUSC_dmrs_d4.rds'))
setDT(dmrs4)
dmrs4_BRCA_gr <- dt2gr(dmrs4[diff.Methy > 0]) # group1=LUSC > group2=BRCA ==> sign. meth. in LUSCC, i.e. BRCA loss of meth
dmrs4_BRCA_gr <- gr.chr(dmrs4_BRCA_gr)
dmrs4_BRCA_gt <- gTrack(dmrs4_BRCA_gr, col='black')
dmrs4_LUSC_gr <- dt2gr(dmrs4[diff.Methy < 0])
dmrs4_LUSC_gr <- gr.chr(dmrs4_LUSC_gr)
dmrs4_LUSC_gt <- gTrack(dmrs4_LUSC_gr, col='black')

# lin.-spec. H3K27ac/ATAC peaks
XX

# ATAC tracks + sign. peaks
IDs <- fread('samplesheet_ATAC_blank.csv',header=T)
IDs_LUSC <- IDs[IDs$tissue=="LUSC",]
IDs_HRpos <- IDs[IDs$tissue=="HR+BRCA",]

# H3K27ac tracks + sign. peaks
IDs <- fread('samplesheet_H3K27ac_blank.csv',header=T)
IDs_LUSC <- IDs[IDs$tissue=="LUSC",]
IDs_HRpos <- IDs[IDs$tissue=="HR+BRCA",]

chip_track_list <- function(FEpath,peakPath) {
  FE <- gTrack(FEpath,bars=T)
  narrowPeaks <- fread(peakPath)
  colnames(narrowPeaks)[1:3] <- c("chr","start","end")
  narrowPeaks_gt <- gTrack(dt2gr(narrowPeaks), height=0.5)
  return(c(FE,narrowPeaks_gt))
}

LUSC_chip_track_list <- lapply(IDs_LUSC$FEpath, IDs_LUSC$peaks, chip_track_list)
names(LUSC_chip_track_list) <- IDs_LUSC$sid
BRCA_chip_track_list <- lapply(IDs_BRCA$FEpath, IDs_BRCA$peaks, chip_track_list)
names(BRCA_chip_track_list) <- IDs_BRCA$sid

# plot
plot(c(gt_ge,
       dmrs4_BRCA_gt,
       #BRCA_meth_track_list,
       walks_consAmp_BRCA_gt,
       sign_LvR_b_gt,
       sign_BvL_gt,
       dmrs4_LUSC_gt,
       #LUSC_meth_track_list,
       walks_consAmp_NSCLC_gt),plot_window)
