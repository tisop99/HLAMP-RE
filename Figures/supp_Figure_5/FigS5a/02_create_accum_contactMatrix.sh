cooler cload pairs -c1 2 -p1 3 -c2 5 -p2 6 \
	../../../common/data/refGenome/T2T_contig_lengths_UCSC.tsv:5000 \
        ${gene}_${tissue}_pooled.allValidPairs \
        ${gene}_${tissue}_pooled_5kb_rawCounts.cool

cooler dump --table pixels --range ${chr} ${gene}_${tissue}_pooled_5kb_rawCounts.cool > ${gene}_${tissue}_pooled_5kb_rawCounts_${chr}.tsv
