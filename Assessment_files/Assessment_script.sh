#!/bin/bash

#Set working directory and make file structure
cd ~
mkdir assessment
cd ./assessment
mkdir data results logs meta other
mkdir -p ./data/trimmed_fastq
mkdir -p ~/assessment/data/reference
mkdir ~/assessment/results/fastqc_trimmed_reads
mkdir ~/assessment/results/fastqc_untrimmed_reads
mkdir ~/assessment/data/aligned_data

#2.1.1 Download dependencies (assumes anaconda/miniconda is already installed)
conda install trimmomatic
conda install fastqc
conda install bwa
conda install samtools
conda install freebayes
conda install picard
conda install bedtools
conda install vcflib

#It is assumed that annovar is already installed in ~/assessment/annovar/
cd ./annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 humandb/     #This one is particularly important, as a later step will be to filter by dbSNP

#Download and unzip snpEff
cd ~
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip

#2.1.2 Download files in to 'data' directory (not neccessary for this script - pipeline assumes user inputs the correct data and path when running the bash script)
#cd ./data
#wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
#wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
#wget  https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
#mv ~/assessment/data/NGS0001.R1.fastq.qz ~/assessment/data/NGS0001.R1.fastq.gz
#mv ~/assessment/data/NGS0001.R2.fastq.qz ~/assessment/data/NGS0002.R1.fastq.gz

#download and move reference genome data to reference directory
cd ./assessment
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv hg19.fa.gz ~/assessment/data/reference

#2.2.1 Quality assessment and read trimming
fastqc -t 4 $1 $2

mv ~/assessment/data/*fastqc* ~/assessment/results/fastqc_untrimmed_reads

trimmomatic PE  \
  -threads 4 \
  -phred33 \
  $1 $2 \
  -baseout ~/assessment/data/trimmed_fastq/trimmed_data \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:5

#2.2.2
fastqc -t 4 ~/assessment/data/trimmed_fastq/trimmed_data_1P \
	~/assessment/data/trimmed_fastq/trimmed_data_2P

mv ~/assessment/data/trimmed_fastq/*fastqc* ~/assessment/results/fastqc_trimmed_reads/

#2.3
#Make Genome Index
bwa index ~/assessment/data/reference/hg19.fa.gz

#Alignment with BWA-mem (to create SAM - not used here because SAM file is too large)
#bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1.111.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001\tDT:2017-02-23\tPU:D1375ACXX.1' -I 250,50 \
#~/assessment/data/reference/hg19.fa.gz ~/assessment/data/trimmed_fastq/trimmed_data_1P \
# ~/assessment/data/trimmed_fastq/trimmed_data_2P > \
# ~/assessment/data/aligned_data/NGS0001.sam

#Alignment with BWA-mem (pipe directly to samtools view to convert to BAM due to disk space limitations)
bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1.111.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001\tDT:2017-02-23\tPU:D1375ACXX.1' -I 212,31 \
~/assessment/data/reference/hg19.fa.gz ~/assessment/data/trimmed_fastq/trimmed_data_1P \
~/assessment/data/trimmed_fastq/trimmed_data_2P | samtools view -h -b >  ~/assessment/data/aligned_data/NGS0001.bam

#Convert sam > bam (completed using pipe above)
 #samtools view -h -b NGS0001.sam > NGS0001.bam

#Sort bam file and index the sorted file
cd ~/assessment/data/aligned_data
 samtools sort NGS0001.bam > NGS0001_sorted.bam
 samtools index NGS0001_sorted.bam

 #Duplicate Marking
 picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt

 #Index sorted and marked bam file
 samtools index NGS0001_sorted_marked.bam

 #Filter BAM based on mapping quality and bitwise flags w/ samtools
 samtools view \
 -F 1796  \
 -q 20 \
 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam

 samtools index NGS0001_sorted_filtered.bam

#Post alignment QC
#View the bam file with: samtools view NGS0001_sorted_filtered.bam | less -S
samtools flagstat NGS0001_sorted_filtered.bam > flagstats_out.txt
samtools idxstats NGS0001_sorted_filtered.bam > alignstats_per_chr.txt
picard CollectInsertSizeMetrics I= NGS0001_sorted_filtered.bam O= insert_size_metrics.txt H= insert_size_histogram.pdf M= 0.05
bedtools coverage \
-a NGS0001_sorted_filtered.bam \
-b ~/assessment/data/annotation.bed  > coverage.txt

#2.4 Freebayes variant calling
zcat ~/assessment/data/reference/hg19.fa.gz > ~/assessment/data/reference/hg19.fa
samtools faidx ~/assessment/data/reference/hg19.fa
freebayes \
--bam ~/assessment/data/aligned_data/NGS0001_sorted_filtered.bam \
--fasta-reference ~/assessment/data/reference/hg19.fa \
--vcf ~/assessment/results/NGS0001.vcf
#--targets ~/assessment/data/annotation.bed  Alternative to limit analysis to only positions in bed file

#Compress and index VCF
bgzip ~/assessment/results/NGS0001.vcf
tabix -p vcf ~/assessment/results/NGS0001.vcf.gz

#Filter VCF
vcffilter \
-f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
	~/assessment/results/NGS0001.vcf.gz > ~/assessment/results/NGS0001_filtered.vcf

#Restrict to only regions in bed file with bedtools intersect
bedtools intersect \
-header -wa -a ~/assessment/results/NGS0001_filtered.vcf \
-b ../annotation.bed > ~/assessment/results/NGS0001_filtered_anno.vcf

#Compress and index new VCF, as before
bgzip ~/assessment/results/NGS0001_filtered_anno.vcf
tabix -p vcf ~/assessment/results/NGS0001_filtered_anno.vcf.gz

#2.5
#VCF - to - annovar input
~/assessment/annovar/convert2annovar.pl \
-format vcf4 ~/assessment/results/NGS0001_filtered_anno.vcf.gz > \
~/assessment/results/NGS0001_filtered_anno.avinput

#Annovar annotation
~/assessment/annovar/table_annovar.pl  \
~/assessment/results/NGS0001_filtered_anno.avinput assessment/annovar/humandb/ \
-buildver hg19 \
-out ~/assessment/results/NGS0001_filtered_anno \
-remove \
-protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro,avsnp150  \
-operation g,g,f,f,f,f \
-otherinfo \
-nastring . \
-csvout

#snpEff Annotation
java -Xmx8g -jar ~/snpEff/snpEff.jar -v -stats snpeff_ann_stats.html hg19 ~/assessment/results/NGS0001_filtered_anno.vcf.gz >  ~/assessment/results/NGS0001_snpeff.ann.vcf
mv ./snpeff_ann_stats* ~/assessment/results/


#java -jar SnpSift.jar annotate -dbsnp ~/assessment/results/NGS0001_snpeff.ann.vcf > ~/assessment/results/NGS0001_filtered_anno.dbsnp.vcf.gz
#java -Xmx8g -jar snpEff.jar hg19 ~/assessment/results/NGS0001_filtered_anno.vcf.gz > NGS0001_snpeff.vcf
