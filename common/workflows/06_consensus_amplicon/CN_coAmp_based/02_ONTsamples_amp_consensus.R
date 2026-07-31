library(gGnome)
library(data.table)

case <- commandArgs(trailingOnly = TRUE)
xG <- case[1] # CCND1, NSD3, EGFR
tt <- case[2] # tissue type (LUSC / BRCA)

IDs_all <- fread(paste0(xG,'_amp_ONTsamples_annotated.csv'),header=T)
IDs <- IDs_all[IDs_all$tissue == tt]
gr_list <- list()

paths <- fread('samplesheet_purity_blank.csv')

for (i in 1:nrow(IDs)) {
  xid <- IDs$sample[i]
  jabba_dir <- paths$jabbaDir[i]

  path_pattern <- paste0(jabba_dir,"/jabba.seg")
  file_path <- Sys.glob(path_pattern)
  if (length(file_path) != 1) { stop("Expected exactly one file, found ", length(file_path)) }
  tmp_dt <- fread(file_path)
  tmp_dt <- tmp_dt[, .(chrom,start,end,cn)]
  tmp_dt <- dt2gr(tmp_dt)
  tmp_dt <- gr.chr(tmp_dt)
  tmp_dt <- tmp_dt[tmp_dt$cn >= 10]   # all highly co-amplified regions
  gr_list <- append(gr_list, list(tmp_dt))
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

saveRDS(collapsed,paste0('ONTsamples_uniStrand_ampSum_CNbased_',xG,'_',tt,'.rds'))
