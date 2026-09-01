#!/bin/bash
#
#SBATCH --job-name=Fig.5A_20260807
#SBATCH --output=/data/cephfs-1/work/projects/dubois-lrwgs/ying/Fig.5A/Fig_5A_CCND1_TF.out
#SBATCH --error=/data/cephfs-1/work/projects/dubois-lrwgs/ying/Fig.5A/Fig_5A_CCND1_TF.err
#SBATCH --ntasks=12
#SBATCH --time=0-1:00:00
#SBATCH --mem=48G

mkdir -p /data/cephfs-1/work/projects/dubois-lrwgs/ying/Fig.5A

# container + resources
# /runins/  = /data/cephfs-2/unmirrored/groups/dubois/ (CGC, gencode, chain, TF ChIP-seq, HMF&PCAWG ampdists)
# /inputdata/ = /data/cephfs-1/work/projects/dubois-lrwgs/frank/ (NSD3 diffmethyl, Fisher, v4C, pdf output)
apptainer exec --bind /data/cephfs-2/unmirrored/groups/dubois/:/runins/ --bind /data/cephfs-1/work/projects/dubois-lrwgs/:/inputdata/ \
    /data/cephfs-2/unmirrored/groups/dubois/jabba_gxg.sif bash -c "

Rscript /inputdata/ying/Fig.5A/20260810_Fig5a_1.r
"
