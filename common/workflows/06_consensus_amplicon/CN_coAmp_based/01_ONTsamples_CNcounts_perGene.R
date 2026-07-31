library(gGnome)
library(data.table)
library(rtracklayer)

CGC_gr_T2T <- readRDS('../../../data/refGenome/cgc_T2T_gr_UCSC.rds')
NSD3 <- CGC_gr_T2T[mcols(CGC_gr_T2T)$Gene_Symbol == 'NSD3']
CCND1 <- CGC_gr_T2T[mcols(CGC_gr_T2T)$Gene_Symbol == 'CCND1']
EGFR <- CGC_gr_T2T[mcols(CGC_gr_T2T)$Gene_Symbol == 'EGFR']
targets <- c(CCND1,EGFR,NSD3)

IDs <- fread('samplesheet_purity_blank.csv')
gr_list <- list()

# get CN per sample
for (i in 1:nrow(IDs)) {
  xid <- IDs$sid[i]
  jabba_dir <- IDs$jabbaDir[i]

  path_pattern <- paste0(jabba_dir,"/jabba.seg")
  file_path <- Sys.glob(path_pattern)
  if (length(file_path) != 1) { stop("Expected exactly one file, found ", length(file_path)) }
  tmp_dt <- fread(file_path)
  tmp_dt <- tmp_dt[, .(chrom,start,end,cn)]
  tmp_dt <- dt2gr(tmp_dt)
  tmp_dt <- gr.chr(tmp_dt)
  subset_gr <- subsetByOverlaps(tmp_dt, targets)
  gr_list[[xid]] <- subset_gr
}

all_dt <- rbindlist(
  lapply(names(gr_list), function(sid) {
    dt <- as.data.table(gr_list[[sid]])
    dt[, chr := as.character(seqnames)]
    dt[, seqnames := NULL][, start := NULL][, end := NULL]  # optional
    dt[, sample := sid]
}
))
wide_dt <- dcast(all_dt, chr ~ sample, value.var = "cn", fun.aggregate = mean)
setnames(wide_dt, "chr", "gene")
chr_map <- c(chr11 = "CCND1", chr7 = "EGFR", chr8 = "NSD3")
wide_dt[, gene := chr_map[gene]]
fwrite(wide_dt, file="ONT_perSample_CN_table_CCND1_EGFR_NSD3.csv")

# get all amp. samples per gene
df_all <- read.csv("ONT_perSample_CN_table_CCND1_EGFR_NSD3.csv", row.names=1,check.names=F)

for (gene in rownames(df_all)) {
  cn_values <- df_all[gene, ]
  selected_samples <- names(cn_values)[cn_values >= 10]
  out <- data.frame(sample=selected_samples)
  write.csv(out,file=paste0(gene,"_amp_ONTsamples.csv"),row.names=F)
}

# annotate tissue type
meta <- fread("ONT_samples_meta.csv",header=T)
for (gene in c("CCND1","NSD3","EGFR")) {
  ampList <- fread(paste0(gene,"_amp_ONTsamples.csv"))
  ampList_annot <- meta[ampList, on = .(sample = sample), .(sample,tissue)]
  write.csv(ampList_annot,file=paste0(gene,"_amp_ONTsamples_annotated.csv"),row.names=F)
}
