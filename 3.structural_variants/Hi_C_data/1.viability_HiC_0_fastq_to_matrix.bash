# make a usable HiC matrix for each subspecies


##STEPS:
#1	Align sequences with bwa
#2	Create cut sites file with hicFindRestSite
#3	Make hic matric with hicBuildMatrix



REF="/scratch/user/sblain/get_AIMs/define_AIMs/inland_ref/inland_genome.fasta"
cd /scratch/user/sblain/misc_tasks/HiC/process_HiCExplorer/ 

#Step one: align data

module purge
module load GCC/13.2.0 SAMtools/1.19.2

fastqDir="/scratch/user/sblain/misc_tasks/HiC/fastqs"

bwa mem -A1 -B4  -E50 -L0 $REF \
    $fastqDir/DTG-HiC-872_R1_001.fastq.gz 2>>CoastalHiC_R1.log | samtools view -Shb - > CoastalHiC_R1.bam

bwa mem -A1 -B4  -E50 -L0 $REF \
    $fastqDir/DTG-HiC-872_R2_001.fastq.gz 2>>CoastalHiC_R2.log | samtools view -Shb - > CoastalHiC_R2.bam

bwa mem -A1 -B4  -E50 -L0 $REF \
    $fastqDir/bCatUst1_ARI1_005_R1.fastq.gz 2>>InlandHiC_R1.log | samtools view -Shb - > InlandHiC_R1.bam
	
bwa mem -A1 -B4  -E50 -L0 $REF \
    $fastqDir/bCatUst1_ARI1_005_R2.fastq.gz 2>>InlandHiC_R2.log | samtools view -Shb - > InlandHiC_R2.bam


samtools sort CoastalHiC_R1.bam -o CoastalHiC_R1.sorted.bam
samtools sort CoastalHiC_R2.bam -o CoastalHiC_R2.sorted.bam
samtools index CoastalHiC_R1.sorted.bam
samtools index CoastalHiC_R2.sorted.bam

samtools sort InlandHiC_R1.bam -o InlandHiC_R1.sorted.bam
samtools sort InlandHiC_R2.bam -o InlandHiC_R2.sorted.bam
samtools index InlandHiC_R1.sorted.bam
samtools index InlandHiC_R2.sorted.bam

#Step two: find restriction sites

module purge
module load GCC/11.3.0 OpenMPI/4.1.4 HiCExplorer/3.7.2 #load HiCExplorer

##Notes on restriction enzymes:
#For coastal - used 4 cutter DpnII (cut site GATC, overhang GATC)
#For inland - Arima - "restriction enzymes that digest chromatin at ^GATC and G^ANTC" (DpnII and HinfI)
#helpful for finding cut sites and overhangs: https://enzymefinder.neb.com/

hicFindRestSite --fasta $REF --searchPattern GATC -o inlandGenome_GATC_cutSites.bed
#search pattern is regexp - '.' represents any base, not 'N'
hicFindRestSite --fasta $REF --searchPattern GA.TC -o inlandGenome_GANTC_cutSites.bed

#merge bed files for inland - probably inefficient but it works
cp inlandGenome_GATC_cutSites.bed inlandGenome_GATC_GANTC_cutSites.tmp.bed
cat inlandGenome_GANTC_cutSites.bed >> inlandGenome_GATC_GANTC_cutSites.tmp.bed #combine into one file
sort -k1,1 -k2,2n inlandGenome_GATC_GANTC_cutSites.tmp.bed > inlandGenome_GATC_GANTC_cutSites.tmp.sorted.bed #sort file by chromosome and start position
bedtools merge -i inlandGenome_GATC_GANTC_cutSites.tmp.sorted.bed > inlandGenome_GATC_GANTC_cutSites.merged.bed #merge intervals from one file
rm inlandGenome_GATC_GANTC_cutSites.tmp.* #clean up files

############### Step three: build HiC Matrix

module purge
module load GCC/11.3.0 OpenMPI/4.1.4 HiCExplorer/3.7.2 #load HiCExplorer


# hicBuildMatrix --samFiles InlandHiC_R1.sorted.bam InlandHiC_R2.sorted.bam \
	# --binSize [10000 - small bin now - can be merged later] \
	# --restrictionSequence [GATC - restriction site(s) - space delimited] \ 
	# --danglingSequence [GATC A.T - overhangs - in same order as restriction sites] \
	# --restrictionCutFile [sorted bed file listing all cut sites] \
	# --threads 4 \
	# --inputBufferSize 100000 \
	# --outBam InlandHiC.bam \
	# -o InlandHiC_matrix.h5 \
	# --QCfolder ./InlandHiC_QC

hicBuildMatrix --samFiles CoastalHiC_R1.sorted.bam CoastalHiC_R2.sorted.bam \
  --binSize 10000 \
  --restrictionSequence GATC \
  --danglingSequence GATC \
  --restrictionCutFile inlandGenome_GATC_cutSites.bed \
  --threads 4 \
  --outBam CoastalHiC.bam \
  -o CoastalHiC_matrix.h5 \
  --QCfolder ./CoastalHiC_QC

hicBuildMatrix --samFiles InlandHiC_R1.sorted.bam InlandHiC_R2.sorted.bam \
  --binSize 10000 \
  --restrictionSequence GATC GA.TC \
  --danglingSequence GATC A.T \
  --restrictionCutFile inlandGenome_GATC_GANTC_cutSites.merged.bed \
  --threads 4 \
  --outBam InlandHiC.bam \
  -o InlandHiC_matrix.h5 \
  --QCfolder ./InlandHiC_QC

############### Step four: correct hic Matrix

#make diagnostic plot to pick z-score threshold
#remove very high outliers (example, 5) and lots of low outliers (ex -2 to -1)
hicCorrectMatrix diagnostic_plot -m CoastalHiC_matrix.h5 -o CoastalHiC_corrected.png
hicCorrectMatrix diagnostic_plot -m InlandHiC_matrix.h5 -o InlandHiC_corrected.png

#Normalize here if multiple samples need to be normazlied to the same read coverage
hicNormalize -m InlandHiC_matrix.h5 CoastalHiC_matrix.h5 \
 --normalize smallest \
 -o InlandHiC.norm.h5 CoastalHiC.norm.h5

#correct HiC matrix
hicCorrectMatrix correct -m CoastalHiC.norm.h5 --filterThreshold -1.5 5 -o CoastalHiC.norm.corrected.h5
hicCorrectMatrix correct -m InlandHiC.norm.h5 --filterThreshold -1.5 5 -o InlandHiC.norm.corrected.h5
