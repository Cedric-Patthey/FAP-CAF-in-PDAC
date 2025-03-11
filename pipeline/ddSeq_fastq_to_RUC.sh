#!/bin/bash -l

#SBATCH -A NAISS2023-22-1275
#SBATCH -p node -n 20
#SBATCH -t 24:00:00
#SBATCH -J ddSeq_fastq_to_ruc

# run as sbatch ddSeq_fastq_to_RUC.sh <in_R1.fastq.gz> <in_R2.fastq.gz> <rootname> <path/to/mappable/genome> <protcod.bed>

# copies the necessary files to scratch
cp ${1} $SNIC_TMP/R1.fastq.gz &
cp ${2} $SNIC_TMP/R2.fastq.gz &
cp -r ${4} $SNIC_TMP/gendir
wait

cp ${5} $SNIC_TMP/bedfile

cp make_ddSeq_RUC_table.py $SNIC_TMP/make_ddSeq_RUC_table.py
cp Whitelist_ddSeq.txt $SNIC_TMP/Whitelist_ddSeq.txt

# takes note of the working directory
pwd=$(pwd)

# moves to scratch
cd $SNIC_TMP

# loads modules
module load bioinfo-tools
module load python
module load star/2.7.11a
module load samtools/1.20
module load BEDTools/2.31.1
module load cd-hit/4.8.1

# readID is renamed in the R2 fastq file
zcat R2.fastq.gz | awk 'NR%4==1{print "@R"(NR+3)/4} NR%4!=1{print $1}' > rn_R2.fq

# maps with star and annotates with bedtools
STAR --runThreadN 18 --readFilesCommand cat --outFileNamePrefix R2_ --readFilesIn rn_R2.fq --genomeDir gendir --outSAMunmapped None --outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5

samtools view -h -b -F4 -F256 R2_Aligned.out.sam | samtools sort -O bam > mapped.bam
echo "mapping completed"

bedtools intersect -a mapped.bam -b bedfile -wo -split -s -f 0.25 -bed | awk '{print $4"\t"$16"\t"$19}' > read_feature.txt
echo "feature annotation completed"

# runs python script to generate the raw umi count
python make_ddSeq_RUC_table.py $3 R1.fastq.gz read_feature.txt

gzip ${3}_RUC_long.txt

# copies output file
cp -r ${3}_RUC_long.txt.gz $pwd
echo "Raw UMI count pipeline completed"