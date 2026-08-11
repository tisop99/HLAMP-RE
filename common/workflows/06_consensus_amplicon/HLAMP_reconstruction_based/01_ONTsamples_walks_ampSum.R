library(gGnome)
library(data.table)

case <- commandArgs(trailingOnly = TRUE)
xG <- case[1]  # CCND1, NSD3, EGFR
tissue <- case[2] # BRCA, NSCLC

CGC <- readRDS('../../../data/refGenome/cgc_T2T_gr_UCSC.rds')
gene <- CGC[mcols(CGC)$Gene_Symbol == xG]
target <- GRanges(seqnames=as.character(seqnames(gene[1])),ranges=IRanges(start=min(start(gene)),end=max(end(gene))))

IDs <- fread('samplesheet_walks_blank.csv')
gr_list <- list()

for (i in 1:nrow(IDs)) {
  xid <- IDs$sid[i]
  walkPath <- IDs$walksPath[i]
  
  walkX <- readRDS(walkPath)
  walk_gr <- unlist(walkX$grl)
  walk_gr <- gr.chr(walk_gr)
  strand(walk_gr) <- "*"

  hits <- findOverlaps(target,walk_gr,type="within")
  if (length(hits)>0 & all(walk_gr[subjectHits(hits)]$cn >= 10)) {
    gr_list <- append(gr_list, list(walk_gr))
  }
}

gr_list <- lapply(gr_list, function(gr) {
  seqlengths(gr) <- NA
  gr
})
grl <- GRangesList(gr_list)
all_gr <- unlist(grl, use.names = FALSE)
gr_by_chr <- split(all_gr, seqnames(all_gr))

atomic_list <- lapply(gr_by_chr, function(gr_chr) {
  breaks <- sort(unique(c(start(gr_chr), end(gr_chr) +1)))

  GRanges(seqnames=unique(seqnames(gr_chr)),
	  ranges=IRanges(start=head(breaks,-1),
                         end=tail(breaks,-1)-1
			 ))
})
# recombine into a single granges object
atomic_ranges <- unlist(GRangesList(atomic_list), use.names = FALSE)


sample_matrix <- sapply(gr_list, function(gr) {
  countOverlaps(atomic_ranges, gr) > 0			 
})
support <- rowSums(sample_matrix)
mcols(atomic_ranges)$sample_support <- support
#

### collapse adjacent ranges that have the same sample_support (currently distinct ranges since different samples contributed to sample_support)
# prepare the is_new_group() function
seqn   <- as.vector(seqnames(atomic_ranges))
starts <- start(atomic_ranges)
ends   <- end(atomic_ranges)
supp   <- mcols(atomic_ranges)$sample_support
# adjacent ranges are assigned to the same group if: seqnames equal & sample_support equal & the 2 ranges are adjacent
is_new_group <- c(
  TRUE,
  seqn[-1] != seqn[-length(seqn)] |
  supp[-1] != supp[-length(supp)] |
  starts[-1] != ends[-length(ends)] + 1
)
group_id <- cumsum(is_new_group)
gr_split <- split(atomic_ranges, group_id)
# now collapse
gr_split <- split(atomic_ranges, group_id)
collapsed_list <- lapply(gr_split, function(gr) {
  out <- reduce(gr)
  mcols(out)$sample_support <- unique(mcols(gr)$sample_support)
  out
})
collapsed <- unlist(GRangesList(collapsed_list), use.names = FALSE)

saveRDS(collapsed,paste0('ONTsamples_uniStrand_HLAMPwalks_ampSum_',xG,'_',tissue,'.rds'))

