#!/usr/bin/env bash


PATH_FILE="data/mappability/"

cd $PATH_FILE

# Merge unmappable regions
# -----------------------------------
zcat Duke_DAC_exclude_bed.gz hg19_low_complexity_regions.gz | cut -f1,2,3 | sort -k1,1 -k2,2n | mergeBed -i stdin | gzip > merged_unmappable_regions.gz


# Remove unmappable regions from CRG high mappable regions
subtractBed -a wgEncodeCrgMapabilityAlign36mer_score1.bed.gz -b merged_unmappable_regions.gz | gzip > hg19.mappable.gz


# Download the chromosome sizes
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes


# get all genome that is unmappable
awk '{OFS="\t"}{print $1, 0, $2}' hg19.chrom.sizes > hg19.chrom.sizes.bed
subtractBed -a  hg19.chrom.sizes.bed  -b hg19.mappable.gz | gzip > hg19.notmappable.gz

# get all genome without the "chr" prefix
zcat hg19.mappable.gz | sed 's/chr//g' | gzip > hg19.mappable.nochr.gz
