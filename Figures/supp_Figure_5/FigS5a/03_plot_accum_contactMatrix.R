library(GxG)
library(gGnome)
library(rtracklayer)

case <- commandArgs(trailingOnly = TRUE)
xG <- case[1]  # CCND1, NSD3, EGFR
chrom <- case[2]  # chr11, chr8, chr7

binsize=5e3
xT_bins <- fread('source_data/normalized_abs.bed')
colnames(xT_bins) <- c('chr', 'start', 'end')
xT_bins[, start := start +1, by = start]

# load accumulated matrices
chr_B <- fread(paste0(xG,'_BRCA_pooled_5kb_rawCounts_',chrom,'.tsv'))
colnames(chr_B) <-  c("i","j","value")
chr_1stBinIndex <- min(chr_B$i)-1
chr_B$i <- chr_B$i - chr_1stBinIndex
chr_B$j <- chr_B$j - chr_1stBinIndex
interactions_dts_chr_B <- chr_B
gmat_BRCA <- gM(gr=dt2gr(xT_bins[chr==chrom]),dat=interactions_dts_chr_B,fill=0,agg.fun=min)

chr_B <- fread(paste0(xG,'_LUSC_pooled_5kb_rawCounts_',chrom,'.tsv'))
colnames(chr_L) <-  c("i","j","value")
chr_1stBinIndex <- min(chr_L$i)-1
chr_L$i <- chr_L$i - chr_1stBinIndex ##366432
chr_L$j <- chr_L$j - chr_1stBinIndex ##366432
interactions_dts_chr_L <- chr_L
gmat_LUSC <- gM(gr=dt2gr(xT_bins[chr==chrom]),dat=interactions_dts_chr_L,fill=0,agg.fun=min)

# gene track
gencode_T2T_gt <- readRDS('../../../common/data/refGenome/T2T_gencode_gt_UCSC.rds')

window <- readRDS(paste0('source_data/',xG,'_plottingWindow.rds'))

pdf(file=paste0(xG,'_hic_accum.pdf'), paper = 'a4', height = 24)
  plot(c(gencode_T2T_gt,
	 gmat_BRCA$gtrack(height=300,clim=c(0,70))),
         window, cex.label=0.1,xaxis.interval=0.5e6, xaxis.unit=1e6, xaxis.suffix="Mb",yaxis.cex=0.6)
  plot(c(gencode_T2T_gt,
	 gmat_LUSC$gtrack(height=300,clim=c(0,70))),
         window, cex.label=0.1,xaxis.interval=0.5e6, xaxis.unit=1e6, xaxis.suffix="Mb",yaxis.cex=0.6)
dev.off()                                                                                                                                                                       

