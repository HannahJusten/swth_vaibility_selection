module purge
module load GCC/11.3.0 OpenMPI/4.1.4 HiCExplorer/3.7.2 #load HiCExplorer

cd /scratch/user/sblain/HiC/compare_HiCexplorer

hicMergeMatrixBins -m ../process_HiCExplorer/CoastalHiC.norm.corrected.h5 -o CoastalHiC.norm.corrected.1Mb.h5 -nb 100 #make 1000000bp bins
hicMergeMatrixBins -m ../process_HiCExplorer/InlandHiC.norm.corrected.h5 -o InlandHiC.norm.corrected.1Mb.h5 -nb 100 #make 1000000bp bins

############### run matrix comparison

hicCompareMatrices --matrices CoastalHiC.norm.corrected.1Mb.h5 InlandHiC.norm.corrected.1Mb.h5 \
 --operation log2ratio \
 -o compareCoastalInland.1Mb.h5

############### make plots to compare coastal vs inland interactions

# chromosome 5
hicPlotMatrix -m compareCoastalInland.1Mb.h5 \
 --clearMaskedBins \
 --region scaffold_4_arrow_ctg1 \
 -o compareCoastalInland.1Mb.chrom5.png

# chromosome 3
hicPlotMatrix -m compareCoastalInland.1Mb.h5 \
 --clearMaskedBins \
 --region super_scaffold_3 \
 -o compareCoastalInland.1Mb.chrom3.png
 
# chromosome 3 inversion
hicPlotMatrix -m compareCoastalInland.1Mb.h5 \
 --clearMaskedBins \
 --region super_scaffold_3:5000000-18000000 \
 -o compareCoastalInland.1Mb.chrom3_inv8-13.png
 
# chromosome 5
hicPlotMatrix -m compareCoastalInland.1Mb.h5 \
 --clearMaskedBins \
 --region scaffold_4_arrow_ctg1:10000000-60000000 \
 -o compareCoastalInland.1Mb.chrom5_10-60Mb.png

############### make plots to check location of centromere

# Chromosome 5 - two methods
hicPlotMatrix -m InlandHiC.norm.corrected.1Mb.h5 \
 --clearMaskedBins \
 --region scaffold_4_arrow_ctg1 \
 --vMax 3000 \
 --log1p True \
 -o HiCInland.1Mb.chrom5.png
 
hicPlotMatrix -m InlandHiC.norm.corrected.1Mb.h5 \
 --clearMaskedBins \
 --region scaffold_4_arrow_ctg1 \
 --log1p \
 -o HiCInland.1Mb.log.chrom5.png
 
# Chromosome 1 - two methods 
hicPlotMatrix -m InlandHiC.norm.corrected.1Mb.h5 \
 --clearMaskedBins \
 --region super_scaffold_1 \
 --vMax 3000 \
 -o HiCInland.1Mb.chrom1.png

hicPlotMatrix -m InlandHiC.norm.corrected.1Mb.h5 \
 --clearMaskedBins \
 --region super_scaffold_1 \
 --log1p \
 -o HiCInland.1Mb.log.chrom1.png
