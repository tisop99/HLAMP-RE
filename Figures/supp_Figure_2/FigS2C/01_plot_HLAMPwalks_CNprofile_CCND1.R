library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(gGnome)

# load profiles
cn_consAmp_BRCA <- readRDS('source_data/cnBased_consensus_ampSum_extendedCohort_CCND1_BRCA_10kb_binned_interpolated.rds')
cn_consAmp_NSCLC <- readRDS('source_data/cnBased_consensus_ampSum_extendedCohort_CCND1_NSCLC_10kb_binned_interpolated.rds')
walks_consAmp_BRCA <- readRDS('source_data/HLAMPwalks_consensus_ampSum_CCND1_BRCA_10kb_binned_interpolated.rds')
walks_consAmp_NSCLC <- readRDS('source_data/HLAMPwalks_consensus_ampSum_CCND1_NSCLC_10kb_binned_interpolated.rds')
CN_B <- as.data.frame(cn_consAmp_BRCA)
CN_L <- as.data.frame(cn_consAmp_NSCLC)
walks_B <- as.data.frame(walks_consAmp_BRCA)
walks_L <- as.data.frame(walks_consAmp_NSCLC)

# add offset to chromosome start (1 + offset)
chr_sizes <- fread('../../../common/data/refGenome/T2T_contig_lengths_UCSC.tsv')
colnames(chr_sizes) <- c("seqnames","chr_len")
chr_sizes <- chr_sizes %>%
  mutate(offset = lag(cumsum(as.numeric(chr_len)), default = 0))

CN_B_plot <- CN_B %>%
  left_join(chr_sizes, by="seqnames") %>%
  mutate(genome_pos=start+offset)

CN_L_plot <- CN_L %>%
  left_join(chr_sizes, by="seqnames") %>%
  mutate(genome_pos=start+offset)

walks_B_plot <- walks_B %>%
  left_join(chr_sizes, by="seqnames") %>%
  mutate(genome_pos=start+offset)

walks_L_plot <- walks_L %>%
  left_join(chr_sizes, by="seqnames") %>%
  mutate(genome_pos=start+offset)

CN_B_plot$seqnames <- factor(CN_B_plot$seqnames, levels=paste0("chr", c(1:22, "X", "M")))
CN_L_plot$seqnames <- factor(CN_L_plot$seqnames, levels=paste0("chr", c(1:22, "X", "M")))
walks_B_plot$seqnames <- factor(walks_B_plot$seqnames, levels=paste0("chr", c(1:22, "X", "M")))
walks_L_plot$seqnames <- factor(walks_L_plot$seqnames, levels=paste0("chr", c(1:22, "X", "M")))

chr_bounds <- CN_B_plot %>%
  distinct(seqnames,chr_len) %>%
  mutate(xmin=0,xmax=chr_len)

# plot
make_plot <- function(data, color) {
  ggplot(data %>% filter(AmpFraction_smoothed > 0)) +
    geom_rect(aes(xmin=start,xmax=end,ymin = 0,ymax=AmpFraction_smoothed),color=color) +
    geom_rect(data=chr_bounds,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
	      inherit.aes=F,fill=NA,color=NA) +
    scale_x_continuous(
      labels = function(x) paste0(round(x/1e6), "Mb"),
      breaks = function(x) c(min(x),max(x)),
      expand = c(0,0)
    ) +
    facet_grid(~seqnames, scales="free_x", space="free_x") +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    labs(x = "Chromosome", y = "Amplification fraction (% samples)", fill = NULL) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(angle=45,hjust=1,size=8),
      axis.title=element_text(face="bold"),
      panel.spacing.x = unit(0.5, "lines"),
      strip.background = element_blank(),
    )
}

t1 <- make_plot(CN_B_plot, "brown3")
t2 <- make_plot(CN_L_plot, "royalblue3")
t3 <- make_plot(walks_B_plot, "brown3")
t4 <- make_plot(walks_L_plot, "royalblue3")

t3 / t4 / t1 / t2
