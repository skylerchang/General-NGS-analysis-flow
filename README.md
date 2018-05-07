# General-NGS-analysis-flow

Part of NGS pipeline analysis example VCF analysis project

# Login to sharcnet and download the genome and VCF files to my directory#

# Gunzip two files in variants directory #

gunzip Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz
gunzip Drosophila_melanogaster.vcf.gz

# Download the NCBI SRA Drosophila melanogaster data set SRR3931619 using fastq-dump with SRAtoolkit #

fastq-dump --split-files --gzip -O reads/ SRR3931619

# Login to Development node to run the following java steps #

ssh orc-dev2

# Move into my Binf6210’s dna directory#

cd /work/hchang02/Binf6210/dna

# Run FastQC program to analyze the quality of downloaded data SRR3931619#

fastqc -noextract -o qc reads/SRR3931619_1.fastq.gz reads/SRR3931619_2.fastq.gz

#Trim sequence reads to improve the quality with Trim Galore program #

trim_galore -o qc --paired --illumina reads/SRR3931619_1.fastq.gz reads/SRR3931619_2.fastq.gz

# Run FastQC on the trimmed data again #

fastqc -noextract -o qc qc/SRR3931619_1_val_1.fq.gz qc/SRR3931619_2_val_2.fq.gz

samtools faidx ref/Drosophila_melanogaster.BDGP6.dna.toplevel.fa 

# Create a BWA index of the reference downloaded Drosophila_melanogaster fasta #

bwa index ref/Drosophila_melanogaster.BDGP6.dna.toplevel.fa



# Align the reads to the reference Drosophila_melanogaster fasta and save the result as DM.sam in aligned directory#

bwa mem -R '@RG\tID:1\tSM:ecoli\tLB:lib1\tPL:illumina\tPU:unit1' ref/Drosophila_melanogaster.BDGP6.dna.toplevel.fa reads/SRR3931619_1_val_1.fq.gz reads/SRR3931619_2_val_2.fq.gz > aligned/DM.sam

# Transfer the sam format file to bam format file with samtools program as DM_sorted.bam #

samtools sort -O bam -o aligned/DM_sorted.bam -T aligned/DM_temp aligned/DM.sam

# Create a samtools index of the new bam file #

samtools index aligned/DM_sorted.bam

# Create another format of the reference fasta for subsequent software #

samtools dict ref/Drosophila_melanogaster.BDGP6.dna.toplevel.fa > ref/Drosophila_melanogaster.BDGP6.dna.toplevel.dict

# Copy and change the name of new format fasta name as a end with .fa.dict #

cp ref/Drosophila_melanogaster.BDGP6.dna.toplevel.dict  ref/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.dict

# Use Genome Analysis Toolkit (GATK) to identify potential variant sites in need of realignment #

java -Xmx2g -jar $GATK -T RealignerTargetCreator \
-R ref/Drosophila_melanogaster.BDGP6.dna.toplevel.fa \
-I aligned/DM_sorted.bam \
-o aligned/DM_sorted.intervals

# Use GATK realignment to improve indels by having fewer mismatches and identifying the consensus indel from all reads # 

java -Xmx4g -jar $GATK -T IndelRealigner \
-R ref/Drosophila_melanogaster.BDGP6.dna.toplevel.fa \
-I aligned/DM_sorted.bam \
-targetIntervals aligned/DM_sorted.intervals \
-o aligned/DM_realigned.bam

#  Use Drosophila_melanogaster.vcf file as a known true variants as reference in BQSR model calibration #

java -Xmx4g -jar $GATK -T BaseRecalibrator \
-R ref/Drosophila_melanogaster.BDGP6.dna.toplevel.fa \
-knownSites ref/Drosophila_melanogaster.vcf \
-I aligned/DM_realigned.bam \
-o aligned/DM_recal.table

# Rewrite the quality scores using the data in the DM_recal.table file into a new BAM file.This step uses the recalibration table data produced by BaseRecalibration to recalibrate the quality scores in input.bam, and writing out a new BAM file DM_recal.bam with recalibrated QUAL field values #

java -Xmx2g -jar $GATK -T PrintReads \
-R ref/Drosophila_melanogaster.BDGP6.dna.toplevel.fa \
-I aligned/DM_realigned.bam \
--BQSR aligned/DM_recal.table \
-o aligned/DM_recal.bam

# Use markdulicates function in Picard software to identify duplicate reads and mask such reads, which can avoiding the same evidenced isn’t used multiple times in calling a variant #

java -jar $PICARD MarkDuplicates \
I=aligned/DM_realigned.bam \
O=aligned/DM_marked.bam \
M=aligned/DM_marked_dups_metrics.txt \
CREATE_INDEX=true

# Use samtools to manipulate new created DM _marked.bam and output a binary variant call format as DM.bcf #

samtools mpileup -go variants/DM.bcf -f ref/Drosophila_melanogaster.BDGP6.dna.toplevel.fa aligned/DM_marked.bam

# Extract variants from the samtools mpileup output DM.bcf with bcftools call #

bcftools call -vm --ploidy 1 -o variants/Drosophila_melanogaster.vcf variants/DM.bcf

# Filter out variants with low quality and select quality over 10 then saved as a new file called DM_filtered.vcf # 

bcftools filter -o variants/DM_filtered.vcf -s LOWQUAL -i'%QUAL>10' variants/Drosophila_melanogaster.vcf

