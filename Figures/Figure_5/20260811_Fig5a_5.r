library(data.table)
library(GenomicRanges)
library(GxG)
library(parallel)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)


## Load core annotation
CGC <- readRDS('/runins/me_proj/data/CGC_hg38.rds')
gencode_hg38 <- readRDS('/runins/me_proj/data/gencode_hg38.rds')
hs1ToHg38_chain <- import.chain('/runins/me_proj/data/20251028_lindiffs/hs1ToHg38.over.chain')

# HMF&PCAWG ampdist
BRCA_xGtad_per_amp_gt <- readRDS(paste0("/inputdata/20260512diffmethyl/20260624_Fig.2C_BRCA_CCND1_AmpFreq_track.rds"))
BRCA_xGtad_per_amp_gt_gr <- unlist(liftOver(BRCA_xGtad_per_amp_gt@data[[1]], hs1ToHg38_chain))
BRCA_xGtad_per_amp_gt <- gTrack(data = BRCA_xGtad_per_amp_gt_gr, y.field = "AmpFreq", name = "%amp", col = "red")
BRCA_xGtad_per_amp_gt$bars <- F

BRCA_CCND1_LvsR_sig_LR_gt <- readRDS(paste0("/inputdata/20260512diffmethyl/20260624_Fig.2C_BRCA_CCND1_LvsR_sig_LR_cons80.rds")) # nolint: line_length_linter.
BRCA_CCND1_LvsR_sig_LR_gt_gr <- unlist(liftOver(BRCA_CCND1_LvsR_sig_LR_gt@data[[1]], hs1ToHg38_chain))
BRCA_CCND1_LvsR_sig_LR_gt <- gTrack(data = BRCA_CCND1_LvsR_sig_LR_gt_gr, name = "nlgp", col = "red", y0 = 0, y1 = 1)

BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- readRDS(paste0("/inputdata/20260512diffmethyl/20260624_Fig.2C_BRCAvsNSCLC_CCND1_chr11_sig_track_cons80.rds")) # nolint: line_length_linter.
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr <- unlist(liftOver(BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt@data[[1]], hs1ToHg38_chain))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- gTrack(data = BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr, name = "nlgp", col = "black")

NSCLC_xGtad_per_amp_gt <- readRDS(paste0("/inputdata/20260512diffmethyl/20260624_Fig.2C_NSCLC_CCND1_AmpFreq_track.rds"))
NSCLC_xG_amps_grc <- NSCLC_xGtad_per_amp_gt@data[[1]]
NSCLC_xGtad_amps_hg38gr <- unlist(liftOver(NSCLC_xG_amps_grc, hs1ToHg38_chain))
NSCLC_xGtad_per_amp_gt <- gTrack(data = NSCLC_xGtad_amps_hg38gr, y.field = "AmpFreq", name = "%amp", col = "blue", y0 = 0, y1 = 1)
NSCLC_xGtad_per_amp_gt$bars <- F


xG <- "CCND1"
Oncgr <- CGC %Q% (Gene_Symbol == xG)
window_oi <- Oncgr + 12e5


## TF binding data
MCF7_ESR1_hg38 <- gTrack('/runins/me_proj/data/20250630_BRCA_FGFR1_epigenetics/MCF7_ESR1_hg38_FE_ENCFF237WTX.bigWig'
                         , bars = T, y1 = 100, name = 'ESR1')
MCF7_FOXA1_hg38 <- gTrack('/runins/me_proj/data/20250630_BRCA_FGFR1_epigenetics/MCF7_FOXA1_hg38_FE_ENCFF795BHZ.bigWig'
                          , bars = T, y1 = 100, name = 'FOXA1')
MCF7_GATA3_hg38 <- gTrack('/runins/me_proj/data/20250630_BRCA_FGFR1_epigenetics/MCF7_GATA3_hg38_FE_ENCFF384CPN.bigWig',
                          bars = T, y1 = 100, name = 'GATA3')
MCF7_TRPS1_hg38 <- gTrack('/runins/ying/ChIP-seq/TFs/bigwig/MCF7_TRPS1_hg38_FE.bw',
                          bars = T, y1 = 40, name = 'TRPS1')
MCF7_SPDEF_hg38 <- gTrack('/runins/ying/ChIP-seq/TFs/MCF7-SPDEF-ENCFF114YMT-hg38.bigWig',
                          bars = T, y1 = 40, name = 'SPDEF')
# MCF7_FOXM1_hg38 <- gTrack('/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/MEprojekt/data/20250819regels_cons_motifs_vis/MCF7_FOXM1_hg_38_FE_ENCFF230WOU.bigWig'
#                           , bars = T, y1 = 20, name = 'FOXM1')
MCF7_TFAP2C_hg38 <- gTrack('/runins/ying/ChIP-seq/TFs/bigwig/MCF7_TFAP2C_hg38_FE.bw',
                          bars = T, y1 = 40, name = 'TFAP2C')                    
# T47D_FOXA1_hg38_FE <- gTrack('/runins/me_proj/data/20250630_BRCA_FGFR1_epigenetics/T47D_FOXA1_hg38_FE_ENCFF589DAO.bigWig'
#                              , bars = T, y1 = 40, name = 'T47D_FOXA1')
# T47D_TRPS1_hg38_FE <- gTrack('/runins/ying/ChIP-seq/TFs/bigwig/T47D_TRSP1_hg38_FE.bw'
#                              , bars = T, y1 = 40, name = 'T47D_TRPS1')
MCF7_ESR1_hg38_bed <-fread('/runins/ying/ChIP-seq/TFs/MCF7_ESR1_ENCFF138XTJ.bed')
colnames(MCF7_ESR1_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_ESR1_hg38_bed_gt <- gTrack(dt2gr(MCF7_ESR1_hg38_bed), height = 0.5)
MCF7_FOXA1_hg38_bed <- fread('/runins/ying/ChIP-seq/TFs/MCF7_FOXA1_ENCFF112JVK.bed')
colnames(MCF7_FOXA1_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_FOXA1_hg38_bed_gt <- gTrack(dt2gr(MCF7_FOXA1_hg38_bed), height = 0.5)
MCF7_GATA3_hg38_bed <- fread('/runins/ying/ChIP-seq/TFs/MCF7-GATA3-ENCFF217KJR-hg38.bed')
colnames(MCF7_GATA3_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_GATA3_hg38_bed_gt <- gTrack(dt2gr(MCF7_GATA3_hg38_bed), height = 0.5)
MCF7_TRPS1_hg38_bed <- fread('/runins/ying/ChIP-seq/TFs/results/bwa/merged_library/macs3/MCF7_TRPS1_hg38_peaks.narrowPeak')
colnames(MCF7_TRPS1_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_TRPS1_hg38_bed_gt <- gTrack(dt2gr(MCF7_TRPS1_hg38_bed), height = 0.5)
MCF7_TFAP2C_hg38_bed <- fread('/runins/ying/ChIP-seq/TFs/results/bwa/merged_library/macs3/MCF7_TFAP2C_hg38_peaks.narrowPeak')
colnames(MCF7_TFAP2C_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_TFAP2C_hg38_bed_gt <- gTrack(dt2gr(MCF7_TFAP2C_hg38_bed), height = 0.5)
MCF7_SPDEF_hg38_bed <- fread('/runins/ying/ChIP-seq/TFs/MCF7-SPDEF-ENCFF464GFC-hg38.bed')
colnames(MCF7_SPDEF_hg38_bed)[1:3] <- c('chr', 'start', 'end')
MCF7_SPDEF_hg38_bed_gt <- gTrack(dt2gr(MCF7_SPDEF_hg38_bed), height = 0.5)
# T47D_FOXA1_hg38_bed <- fread('/runins/ying/ChIP-seq/TFs/T47D-FOXA1-hg38-ENCFF420MLJ.bed')
# colnames(T47D_FOXA1_hg38_bed)[1:3] <- c('chr', 'start', 'end')
# T47D_FOXA1_hg38_bed_gt <- gTrack(dt2gr(T47D_FOXA1_hg38_bed), height = 0.5)
# T47D_TRPS1_hg38_FE_bed <- fread('/runins/ying/ChIP-seq/TFs/results/bwa/merged_library/macs3/T47D_TRSP1_hg38_peaks.narrowPeak')
# colnames(T47D_TRPS1_hg38_FE_bed)[1:3] <- c('chr', 'start', 'end')
# T47D_TRPS1_hg38_FE_bed_gt <- gTrack(dt2gr(T47D_TRPS1_hg38_FE_bed), height = 0.5)

# diffmethyl -> no access
## Differential methylation (tumor) – FGFR1-centric but kept as in original
# dmls <- fread('/runins/me_proj/data/20250612_FGFR1epigenetics/20250531_FGFR1_TM_dmls.csv')
# dmls[, end := pos + 1]
# dmls_BRCA_gt <- gTrack(dt2gr(dmls[diff > 0.5]))
# rm(dmls)
# dmrs

dmrs_d5 <- fread('/runins/ying/DSS/20260512_CCND1_TM_dmrs_d.csv')
colnames(dmrs_d5)[6:7] <- c("meanMethy_BRCA", "meanMethy_BRCA")
dmrs_d5_gr_hs1 <- dt2gr(dmrs_d5[diff.Methy > 0.6])
seqlevels(dmrs_d5_gr_hs1) <- paste0("chr", seqlevels(dmrs_d5_gr_hs1))
dmrs_d5_gr_hg38 <- liftOver(dmrs_d5_gr_hs1, hs1ToHg38_chain)
dmrs_d5_gr_hg38 <- unlist(dmrs_d5_gr_hg38)
# check
seqlevelsStyle(dmrs_d5_gr_hg38)
head(dmrs_d5_gr_hg38)
# make track
dmrs_d5_BRCAgt <- gTrack(dmrs_d5_gr_hg38)

# CCND1 BRCA amplification frequency + left-vs-right and BRCA-vs-NSCLC significance tracks (hs1 -> hg38)
BRCA_xGtad_per_amp_gt <- readRDS(paste0("/inputdata/20260512diffmethyl/20260624_Fig.2C_BRCA_CCND1_AmpFreq_track.rds"))
BRCA_xGtad_per_amp_gt_gr <- unlist(liftOver(BRCA_xGtad_per_amp_gt@data[[1]], hs1ToHg38_chain))
BRCA_xGtad_per_amp_gt <- gTrack(data = BRCA_xGtad_per_amp_gt_gr, y.field = "AmpFreq", name = "%amp", col = "red")
BRCA_xGtad_per_amp_gt$bars <- F

BRCA_CCND1_LvsR_sig_LR_gt <- readRDS(paste0("/inputdata/20260512diffmethyl/20260624_Fig.2C_BRCA_CCND1_LvsR_sig_LR_cons80.rds"))
BRCA_CCND1_LvsR_sig_LR_gt_gr <- unlist(liftOver(BRCA_CCND1_LvsR_sig_LR_gt@data[[1]], hs1ToHg38_chain))
BRCA_CCND1_LvsR_sig_LR_gt <- gTrack(data = BRCA_CCND1_LvsR_sig_LR_gt_gr, name = "nlgp", col = "red", y0 = 0, y1 = 1, height = 0.5)

BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- readRDS(paste0("/inputdata/20260512diffmethyl/20260624_Fig.2C_BRCAvsNSCLC_CCND1_chr11_sig_track_cons80.rds"))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr <- unlist(liftOver(BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt@data[[1]], hs1ToHg38_chain))
BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt <- gTrack(data = BRCAvsNSCLC_CCND1_chr11_sig_track_cons80_gt_gr, name = "nlgp", col = "black")

# CCND1 BRCA virtual-4C combined (Fisher) p-values (hs1 -> hg38)
comb_Pvals_v4C_CCND1_BRCA_gr <- readRDS('/inputdata/20260702_v4C_combp_Fisher/20260702_comb_Pvals_v4C_CCND1_BRCA_gr.rds')
comb_Pvals_v4C_CCND1_BRCA_gr_hg38 <- grl.unlist(liftOver(comb_Pvals_v4C_CCND1_BRCA_gr, hs1ToHg38_chain))
comb_Pvals_v4C_CCND1_BRCA_gt <- gTrack(comb_Pvals_v4C_CCND1_BRCA_gr_hg38, y.field = 'nlg_comb_q', bars = T, col = 'red')

# ATAC&H3K27ac Fisher peaks
fishres_A2 <- fread(paste0("/runins/me_proj/data/20250817combined_npWGS_HMF_PCAWG_ampdists/H3K27ac&ATAC_fisherstats/20250817_"
                           ,xG,"_H3K27ac&ATAC_linDiffs_stats.csv"))
fishres_A2[, seqnames := unlist(strsplit(x = peakID, split = '[:]'))[1], by = peakID]
fishres_A2[, start := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[2]), by = peakID]
fishres_A2[, end := as.numeric(unlist(strsplit(x = peakID, split = '[:]'))[3]), by = peakID]
fishres_A2_gr<- dt2gr(fishres_A2)
fishres_A2_sigs <- fishres_A2_gr %Q% (q_val < 0.05)
fishres_A2_sigs_HRBC <- fishres_A2_sigs %Q% (tstdttype == 'BRCA_HR+')
fishres_sigs_BRCA_HR_gt <- gTrack(data = reduce(fishres_A2_sigs_HRBC), name = 'Fish_AT&H3')

## Functional validation datasets
GERCARres <- as.data.table(read.csv2('/runins/me_proj/data/20250421CL_chipseq_macs/TCGA_ATACseq/GERCARres.csv'))
GERCARres_gr <- dt2gr(GERCARres[chromosome == as.character(seqnames(gr.chr(Oncgr)))])
GERCARres_prolifd21_gt <- gTrack(gr.nochr(GERCARres_gr), y.field = 'logFC_prolif_21', bars = TRUE)

GERCARres_gr_sig <- dt2gr(GERCARres[FDR_prolif_21 < 0.1 & chromosome == as.character(seqnames(gr.chr(Oncgr)))])
GERCARres_prolifd21_q0.1_gt <- gTrack(gr.nochr(GERCARres_gr_sig))

SID_MCF7_E2 <- as.data.table(read.csv2('/runins/me_proj/data/20250613_FGFR1_epigenetics/cd-23-1157_supplementary_table_s3_suppst3_MCF7_E2.csv'))
SID_MCF7_E2_gr <- dt2gr(SID_MCF7_E2)
SID_MCF7_E2_prolifd21_gt <- gTrack(SID_MCF7_E2_gr, y.field = 'log2FC_day21_vs_day07', bars = TRUE)
SID_MCF7_E2_sig <- SID_MCF7_E2[FDR_day21_vs_day07 < 0.1]
SID_MCF7_E2_sig_prolifd21_q0.1_gt <- gTrack(dt2gr(SID_MCF7_E2_sig))

SID_T47D_E2 <- as.data.table(read.csv2('/runins/me_proj/data/20250613_FGFR1_epigenetics/cd-23-1157_supplementary_table_s4_suppst4_T47D_E2.csv'))
SID_T47D_E2_gr <- dt2gr(SID_T47D_E2)
SID_T47D_E2_prolifd21_gt <- gTrack(SID_T47D_E2_gr, y.field = 'log2FC_day21_vs_lib', bars = TRUE)
SID_T47D_E2_sig <- SID_T47D_E2[FDR_day21_vs_lib < 0.1]
SID_T47D_E2_sig_prolifd21_q0.1_gt <- gTrack(dt2gr(SID_T47D_E2_sig))


CCND1_TSS <- fread('/runins/tino/merge_vis/CCND1_TSS.csv')
colnames(CCND1_TSS)[4] = 'seqnames'
colnames(CCND1_TSS)[8] = 'start'
CCND1_TSS[, end := start +1]
CCND1_TSSgr <- dt2gr(CCND1_TSS)

window_oiz <- gr.nochr(GRanges("chr11:69000000-70000000") - 2e5)
window_oiz2 <- gr.nochr(GRanges("chr11:69000000-70000000") - 3e5)

pdf(file = paste0('/inputdata/20260803Fig5/20260812_Fig5a_CCND1_TF_BRCA_zo_5e5.pdf')
    ,paper = 'a4', height = 70)
  plot(c(gencode_hg38, GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt,
         dmrs_d5_BRCAgt, comb_Pvals_v4C_CCND1_BRCA_gt,  
         MCF7_TFAP2C_hg38_bed_gt, MCF7_TFAP2C_hg38, MCF7_TRPS1_hg38_bed_gt, MCF7_TRPS1_hg38, MCF7_SPDEF_hg38_bed_gt, MCF7_SPDEF_hg38, MCF7_GATA3_hg38_bed_gt, MCF7_GATA3_hg38, MCF7_FOXA1_hg38_bed_gt, MCF7_FOXA1_hg38, MCF7_ESR1_hg38_bed_gt, MCF7_ESR1_hg38,
         BRCA_CCND1_LvsR_sig_LR_gt, 
         BRCA_xGtad_per_amp_gt)
       , window_oiz, cex.label = 0.1)
  plot(c(gencode_hg38, GERCARres_prolifd21_q0.1_gt, GERCARres_prolifd21_gt,
         dmrs_d5_BRCAgt, comb_Pvals_v4C_CCND1_BRCA_gt,  
         MCF7_TFAP2C_hg38_bed_gt, MCF7_TFAP2C_hg38, MCF7_TRPS1_hg38_bed_gt, MCF7_TRPS1_hg38, MCF7_SPDEF_hg38_bed_gt, MCF7_SPDEF_hg38, MCF7_GATA3_hg38_bed_gt, MCF7_GATA3_hg38, MCF7_FOXA1_hg38_bed_gt, MCF7_FOXA1_hg38, MCF7_ESR1_hg38_bed_gt, MCF7_ESR1_hg38,
         BRCA_CCND1_LvsR_sig_LR_gt, 
         BRCA_xGtad_per_amp_gt)
       , window_oiz2, cex.label = 0.1)
dev.off()
# # 
