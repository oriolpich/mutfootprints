#!/usr/bin/env bash

mkdir data/capecitabine_leishmania

cd data/capecitabine_leishmania/

# Download the L.infantum JPCM5 genome
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-43/fasta/protists_euglenozoa1_collection/leishmania_infantum_jpcm5/dna/Leishmania_infantum_jpcm5.ASM287v2.dna_rm.toplevel.fa.gz
gunzip Leishmania_infantum_jpcm5.ASM287v2.dna_rm.toplevel.fa.gz

# Download from [ENA](https://www.ebi.ac.uk/ena) the fastq sequences to align with bowtie2

# WT
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR174/ERR174230/ERR174230_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR174/ERR174230/ERR174230_2.fastq.gz

# 5-FU treated Resistant
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR248/ERR248703/ERR248703_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR248/ERR248703/ERR248703_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR248/ERR248704/ERR248704_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR248/ERR248704/ERR248704_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR248/ERR248705/ERR248705_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR248/ERR248705/ERR248705_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR248/ERR248706/ERR248706_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR248/ERR248706/ERR248706_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR248/ERR248707/ERR248707_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR248/ERR248707/ERR248707_2.fastq.gz

# Align the fastq files with [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

# Create the index
bowtie2-build Leishmania_infantum_jpcm5.ASM287v2.dna_rm.toplevel.fa Leishmania_infantum_jpcm5

# Run the alignment with 16 threads
# Original
bowtie2 -p 16 -x Leishmania_infantum_jpcm5 -1 ERR174230_1.fastq.gz -2 ERR174230_2.fastq.gz -S ERR174230.sam
# Resistant
bowtie2 -p 16 -x Leishmania_infantum_jpcm5 -1 ERR248703_1.fastq.gz -2 ERR248703_2.fastq.gz -S ERR248703.sam
bowtie2 -p 16 -x Leishmania_infantum_jpcm5 -1 ERR248704_1.fastq.gz -2 ERR248704_2.fastq.gz -S ERR248704.sam
bowtie2 -p 16 -x Leishmania_infantum_jpcm5 -1 ERR248705_1.fastq.gz -2 ERR248705_2.fastq.gz -S ERR248705.sam
bowtie2 -p 16 -x Leishmania_infantum_jpcm5 -1 ERR248706_1.fastq.gz -2 ERR248706_2.fastq.gz -S ERR248706.sam
bowtie2 -p 16 -x Leishmania_infantum_jpcm5 -1 ERR248707_1.fastq.gz -2 ERR248707_2.fastq.gz -S ERR248707.sam

# Sort the alignments with [samtools](http://www.htslib.org/) with 16 threads

# Original
samtools sort -@ 16 -o ERR174230.bam ERR174230.sam
# Resistant
samtools sort -@ 16 -o ERR248703.bam ERR248703.sam
samtools sort -@ 16 -o ERR248704.bam ERR248704.sam
samtools sort -@ 16 -o ERR248705.bam ERR248705.sam
samtools sort -@ 16 -o ERR248706.bam ERR248706.sam
samtools sort -@ 16 -o ERR248707.bam ERR248707.sam

# Create the pileup with [bcftools](http://www.htslib.org/doc/bcftools.html) and do the variant calling as in the original paper
bcftools mpileup -Ou -f Leishmania_infantum_jpcm5.ASM287v2.dna_rm.toplevel.fa ERR174230.bam ERR248703.bam ERR248704.bam ERR248705.bam ERR248706.bam ERR248707.bam | bcftools call -mv -Ob -o leishm_calls_resistant_and_wt.bcf
bcftools view -i '%QUAL>=20' -C 1 leishm_calls_resistant_and_wt.bcf > leishm_calls_resistant_and_wt_vcf.AC_1.vcf

# Go back to original directory and run scripts
cd ../../
python src/capecitabine_leishmania/format_chromosomes_for_bgdata.py data/capecitabine_leishmania/Leishmania_infantum_jpcm5.ASM287v2.dna_rm.toplevel.fa data/capecitabine_leishmania/bgreference/leish/
python src/capecitabine_leishmania/analysis_spectra_treated_WT.py  data/capecitabine_leishmania/leishm_calls_resistant_and_wt_vcf.AC_1.vcf
