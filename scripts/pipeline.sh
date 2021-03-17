#!/bin/bash


trimmomatic PE  \
  -threads 4 \
  -phred33 \
  $1 $2 \
  -baseout ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:5

fastqc -t 4 /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P \
	/home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P

mkdir ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads

mv ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/*fastqc* ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads/

#Make Index (or in this case, just move it to the right place, as we've already done this step!)
mkdir -p ~/ngs_course/dnaseq_pipeline/data/reference

mv ~/ngs_course/dnaseq/data/reference/hg19* ~/ngs_course/dnaseq_pipeline/data/reference/

#bwa index ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz

ls ~/ngs_course/dnaseq_pipeline/data/reference

#Alignment with BWA-mem
mkdir ~/ngs_course/dnaseq_pipeline/data/aligned_data

#This step already run and takes ages, so just move the file!
mv /home/ubuntu/ngs_course/dnaseq/data/aligned_data/WES01_chr22m.sam /home/ubuntu/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m.sam
#bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50 \
#~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/WES01_chr22m_trimmed_R_1P.fastq \
# ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/WES01_chr22m_trimmed_R_2P.fastq > \
# ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m.sam

#Convert sam > bam
 cd ~/ngs_course/dnaseq_pipeline/data/aligned_data
 samtools view -h -b WES01_chr22m.sam > WES01_chr22m.bam

 samtools sort WES01_chr22m.bam > WES01_chr22m_sorted.bam

 samtools index WES01_chr22m_sorted.bam

 #Duplicate Marking
 picard MarkDuplicates I=WES01_chr22m_sorted.bam O=WES01_chr22m_sorted_marked.bam M=marked_dup_metrics.txt

 samtools index WES01_chr22m_sorted_marked.bam

 #Filter BAM based on mapping quality and bitwise flags w/ samtools
 samtools view \
 -F 1796  \
 -q 20 \
 -o WES01_chr22m_sorted_filtered.bam WES01_chr22m_sorted_marked.bam

 samtools index WES01_chr22m_sorted_filtered.bam

 #Post alignment QC
samtools flagstat WES01_chr22m_sorted_filtered.bam > flagstats_out.txt

#samtools view WES01_chr22m_sorted_filtered.bam | less -S
samtools idxstats WES01_chr22m_sorted_filtered.bam > alignstats_per_chr.txt

picard CollectInsertSizeMetrics I= WES01_chr22m_sorted_filtered.bam O= insert_size_metrics.txt H= insert_size_histogram.pdf M= 0.05

bedtools coverage \
-a WES01_chr22m_sorted_filtered.bam \
-b /home/ubuntu/ngs_course/dnaseq_pipeline/data/chr22.genes.hg19.bed > chr_22_coverage.txt

#Freebayes variant calling
zcat ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa

samtools faidx ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa

freebayes \
--bam ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m_sorted_filtered.bam \
--fasta-reference ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa \
--vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf

bgzip ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf

tabix \
-p vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf.gz

#Filter VCF
vcffilter \
-f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
	~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered.vcf

bedtools intersect \
-header -wa -a ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered.vcf \
-b ../chr22.genes.hg19.bed > ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf

bgzip ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf

tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf.gz

#VCF - to - annovar input
~/ngs_course/dnaseq/annovar/convert2annovar.pl \
-format vcf4 ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf.gz > \
~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.avinput

~/ngs_course/dnaseq/annovar/convert2annovar.pl \
~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.avinput humandb/ \
-buildver hg19 \
-out ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22 \
-remove \
-protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro \
-operation g,g,f,f,f \
-otherinfo \
-nastring . \
-csvout
