
mkdir data/megabase_probability
cd data/megabase_probability

# get 1MB windows
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
bedtools makewindows -g hg19.chrom.sizes -w 1000000 |  grep -v "_" | grep -v chrM  > hg19.1Mb.windows.bed
intersectBed -a ../mappability/hg19.mappable.gz  -b hg19.1Mb.windows.bed  -wo | awk '{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $5 "_" $6,$2+1}' | gzip > hg19.mappable.1Mb.windows.bed.extra.gz

