hicproRes <- "../../../common/workflows/08_HiC/results/hicpro/valid_pairs"

# randomly downsample all samples of the same tissue to same pairs count
shuf -n ${minPairs} ${hicproRes}/${sid}.allValidPairs > ${sid}_sub.allValidPairs
# concatenate
cat *_sub.allValidPairs > ${gene}_${tissue}_pooled.allValidPairs
# clean up
rm *_sub.allValidPairs

# downsample the pooled valid pairs file of tissue type with higher pooled raw counts
shuf -n ${minPooledPairs} ${gene}_${tissue}_pooled.allValidPairs > ${gene}_${tissue}_pooled.downSampled.allValidPairs
rm ${gene}_${tissue}_pooled.allValidPairs
mv ${gene}_${tissue}_pooled.downSampled.allValidPairs ${gene}_${tissue}_pooled.allValidPairs
