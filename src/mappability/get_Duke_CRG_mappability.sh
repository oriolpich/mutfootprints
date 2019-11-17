#!/usr/bin/env bash

# Path of the data (change to the desired PATH)

PATH_FILE="data/mappability"
mkdir -p ${PATH_FILE}

cd $PATH_FILE

# UCSC browser blacklisted region
# ----------------------------------------------------------
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz

# merged the two list using mergeBed (from bedtools: https://bedtools.readthedocs.io/en/latest/) with the following
zcat wgEncodeDukeMapabilityRegionsExcludable.bed.gz wgEncodeDacMapabilityConsensusExcludable.bed.gz | sort -k1,1 -k2,2n | mergeBed -i stdin | gzip -c > Duke_DAC_exclude_bed.gz


# CRG mappability data for 36mer
# ---------------------------------------------------------
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig

# conver to bed format with bigWigToBedGraph
bigWigToBedGraph wgEncodeCrgMapabilityAlign36mer.bigWig wgEncodeCrgMapabilityAlign36mer.bed
gzip -9 wgEncodeCrgMapabilityAlign36mer.bed

# get region that are uniquely mapable (with two mismatches allowed)
zcat wgEncodeCrgMapabilityAlign36mer.bed.gz | awk '$4==1' | gzip -c > wgEncodeCrgMapabilityAlign36mer_score1.bed.gz
