####################################################################
# Script run on cluster
####################################################################

#module load GCC/11.3.0 OpenMPI/4.1.4 HiCExplorer/3.7.2

# Chromosome 5
hicPlotMatrix -m compareCoastalInland.1Mb.h5 \
 --clearMaskedBins \
 --region scaffold_4_arrow_ctg1:10000000-60000000 \
 -o compareCoastalInland.1Mb.chrom5_10-60Mb.pdf


 # Chromosome 1
hicPlotMatrix -m compareCoastalInland.1Mb.h5 \
 --clearMaskedBins \
 --region super_scaffold_1:00000000-170000000 \
 -o compareCoastalInland.1Mb.chrom1_0-170Mb.pdf

 ####################################################################
