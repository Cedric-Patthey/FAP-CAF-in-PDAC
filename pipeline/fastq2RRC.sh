#!/bin/bash -l

# Use as: sbatch fastq2RRC.sh <infile.fastq.gz> <rootname> <parameters.csv> <scripts_directory> <threads>
# The output files will be called rootname.bam, e.g s1.bam

#SBATCH -A NAISS2023-22-1275
#SBATCH -p core -n 5
#SBATCH -t 48:00:00
#SBATCH -J fastq2RRC

# parametrers file is parsed
map_gen_dir=$(cat $3 | awk -F "," '$1 ~ /mappable_genome/ {print $2}')
rRNA_fasta_file=$(cat $3 | awk -F "," '$1 ~ /rRNA_fasta/ {print $2}')
annotation_gtf_file=$(cat $3 | awk -F "," '$1 ~ /annotation_gtf/ {print $2}')
annotation_bed_file=$(cat $3 | awk -F "," '$1 ~ /annotation_bed/ {print $2}')
intronic_bed_file=$(cat $3 | awk -F "," '$1 ~ /intronic_bed/ {print $2}')
cDNA_position_key_file=$(cat $3 | awk -F "," '$1 ~ /cDNA_position_key/ {print $2}')
genic_cap=$(cat $3 | awk -F "," '$1 ~ /genic_cap/ {print $2}')
sym_q_cutoff=$(cat $3 | awk -F "," '$1 ~ /symplicity_q_cutoff/ {print $2}')
sym_FC_cutoff=$(cat $3 | awk -F "," '$1 ~ /symplicity_FC_cutoff/ {print $2}')
tp_bias_max_length=$(cat $3 | awk -F "," '$1 ~ /tp_bias_max_length/ {print $2}')
tp_bias_sample_size=$(cat $3 | awk -F "," '$1 ~ /tp_bias_sample_size/ {print $2}')

# copies the necessary files to scratch
cp ${1} $SNIC_TMP/infq.fastq.gz
cp -r $4 $SNIC_TMP/scripts_dir
cp -r $map_gen_dir $SNIC_TMP/gendir
cp $rRNA_fasta_file $SNIC_TMP/rRNA.fas
cp $annotation_gtf_file $SNIC_TMP/annot.gtf
cp $annotation_bed_file $SNIC_TMP/annot.bed
cp $intronic_bed_file $SNIC_TMP/intronic.bed
cp $cDNA_position_key_file $SNIC_TMP/cDNA_pos_key

# takes note of the working directory
pwd=$(pwd)

# moves to scratch
cd $SNIC_TMP

# samples 10k reads
module load bioinfo-tools
module load seqtk
seqtk sample -s42 infq.fastq.gz 10000 > sub10k.fq

# maps the 10k reads
module load star
module load samtools

STAR --runThreadN $5 --outFileNamePrefix sub10k_ --readFilesIn sub10k.fq --genomeDir gendir --outSAMunmapped Within --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66

samtools view -b -F4 -F256 sub10k_Aligned.out.sam | samtools sort -O bam > sub10k.bam

# counts reads mapped to genes for sub10k
module load htseq
htseq-count -f bam -r pos -i gene_id --additional-attr=gene_name -t exon  -m union --stranded=no sub10k.bam  annot.gtf > sub10k_ht_count

# calculates expected genic and capped subsample size
nraw_fastq=$(zcat infq.fastq.gz | awk 'NR%4==1' | wc -l )
ngenic_10k=$(cat sub10k_ht_count | grep -v '__' | awk '{sum+=$3;} END{print sum;}')
capped_subsam_size=$(( genic_cap * 10000 /  ngenic_10k ))

# subsamples with seqtk
if [[ ${nraw_fastq} -gt ${genic_cap} ]]
then
  seqtk sample -s42 infq.fastq.gz ${capped_subsam_size} > capsized.fastq
else
  zcat infq.fastq.gz > capsized.fastq
fi

echo "capping completed"

# adapters are removed with cutadapt
module load cutadapt

# runs cutadapt
cutadapt -g AAGCAGTGGTATCAACGCAGAGTACGGG  -a CCCGTACTCTGCGTTGATACCACTGCTT -g AAGCAGTGGTATCAACGCAGAGTACTTT  -a AAAGTACTCTGCGTTGATACCACTGCTT -g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -n 5 -e 0.2 -O 2 -m 20 -j $5 -o noAd.fastq.gz capsized.fastq

echo "adapter removal completed"

# runs bowtie2 and removes rRNA
module load bowtie2
module load biopython

bowtie2-build rRNA.fas rRNA

bowtie2 -p $5 -x rRNA -U noAd.fastq.gz -S rRNA_map.sam

samtools view -F256 -F4 rRNA_map.sam | cut -f1 | uniq | sort | uniq > rRNA_unwanted_ids

python scripts_dir/get_fastq_not_from_id_list.py noAd.fastq.gz rRNA_unwanted_ids norRNA.fastq

echo "rRNA removal completed"

# Runs STAR
STAR --runThreadN $5 --outFileNamePrefix ${2}_ --readFilesIn norRNA.fastq --genomeDir gendir --outSAMunmapped Within --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66

samtools view -b -F4 -F256 ${2}_Aligned.out.sam | samtools sort -O sam > ${2}.sam
samtools view -b -f4 -F256 ${2}_Aligned.out.sam > ${2}_unmapped.bam

echo "mapping completed"

# duplicate removal
python scripts_dir/symplicity.py $2.sam ${sym_q_cutoff} ${sym_FC_cutoff}

echo "duplicate removal completed"

# spliced reads
samtools view -h ${2}_noDupl.sam | awk '$6 ~ /N/ || $1 ~ /^@/' > ${2}_spliced.sam

echo "spliced reads selection completed"

# 3' bias
python scripts_dir/get_cdna_pos.py ${2}_noDupl.sam cDNA_pos_key annot.bed ${tp_bias_max_length} ${tp_bias_sample_size}
cat sam_cDNA_pos_df | awk -v tagvar=${2} 'NR!=1{print tagvar"\t"$0} NR==1{print "Sample\t"$0}' > ${2}_cDNA_pos_df

echo "Three prime bias analysis completed"

# ht-seq count
htseq-count -f sam -r pos -i gene_id --additional-attr=gene_name -t exon  -m union --stranded=no ${2}_noDupl.sam  annot.gtf > ${2}_ht_count
htseq-count -f sam -r pos -i gene_id --additional-attr=gene_name -t exon  -m union --stranded=no ${2}_spliced.sam  annot.gtf > ${2}_spliced_ht_count

echo "read count completed"


# conting of intronic reads
module load BEDTools/2.31.1
samtools view -h -Sb ${2}_noDupl.sam > ${2}_noDupl.bam
echo Feature$'\t'Feature_size$'\t'Read_Count > ${2}_intronic_count.txt
bedtools coverage -a intronic.bed -b ${2}_noDupl.bam | awk '{print $4"\t"$9"\t"$7}' >> ${2}_intronic_count.txt

echo "intronic count completed"


# Counting of read types
echo Sample$'\t'Raw_fastq$'\t'Capped_fastq$'\t'noAd$'\t'noAd_norRNA$'\t'mapped$'\t'mapped_noDupl$'\t'Genic$'\t'Spliced$'\t'nGenes$'\t'mitoch > ${2}_fastq_bam_genic_count.txt

ncapped_fastq=$(cat capsized.fastq | awk 'NR%4==1' |  wc -l)
nnoad_fastq=$(zcat noAd.fastq.gz | awk 'NR%4==1' |  wc -l)
nnorna_fastq=$(cat norRNA.fastq | awk 'NR%4==1' | wc -l)
nmapped=$(samtools view ${2}.sam | cut -f1 | uniq | sort | uniq | wc -l)
nmapped_noDupl=$(samtools view ${2}_noDupl.sam | cut -f1 | uniq | sort | uniq | wc -l)
ngenic=$(cat ${2}_ht_count | grep -v '__' | awk '{sum+=$3;}END{print sum;}')
nspliced=$(cat ${2}_spliced_ht_count | grep -v '__' | awk '{sum+=$3;}END{print sum;}')
ngenes=$(cat ${2}_ht_count | grep -v '__' | awk '$3>0{print $2}' | sort | uniq | wc -l)
nmitoch=$(cat ${2}_ht_count | grep -v '__' | awk '$2~/mt-/' | awk '{sum+=$3;}END{print sum;}')

echo -e ${2}'\t'$nraw_fastq'\t'$ncapped_fastq'\t'$nnoad_fastq'\t'$nnorna_fastq'\t'$nmapped'\t'$nmapped_noDupl'\t'$ngenic'\t'$nspliced'\t'$ngenes'\t'$nmitoch >> ${2}_fastq_bam_genic_count.txt

echo "read statistics completed"

# Addition of sample tag to RRC file
cat ${2}_ht_count | awk -v varsam=${2} '{print varsam"\t"$0}' | grep -v '__' > ${2}_ht_count.txt

# Saving of Arch directory
mkdir ${2}_Arch
mv norRNA.fastq ${2}_Arch/${2}_Qcleaned_R1.fastq.gz
samtools view -h -Sb ${2}.sam > ${2}.bam
mv ${2}.bam ${2}_Arch
mv ${2}_intronic_count.txt ${2}_Arch
mv ${2}_Log.final.out ${2}_Arch
mv ${2}_SJ.out.tab ${2}_Arch
mv ${2}_unmapped.bam ${2}_Arch

# copies output files back to working directory
cp ${2}_fastq_bam_genic_count.txt $pwd
cp ${2}_ht_count.txt $pwd
cp ${2}_cDNA_pos_df $pwd
cp -r ${2}_Arch $pwd