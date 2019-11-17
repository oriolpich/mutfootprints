#!/usr/bin/env bash

# ANNOVAR folder
ANNOVAR='src/external_software/annovar'
# input file
vcf_file=$1

# location where the file will be stored
outpath=$2
mkdir -p ${outpath}

annot=".annot.txt"
name_file=$(basename ${vcf_file})

# file to remove
intermediate_file=$TMPDIR${name_file}${annot}

# outfile name
gzipext=".annot.gz"
finalfile=${outpath}/${name_file}${gzipext}

# prepare annovar format
perl ${ANNOVAR}/convert2annovar.pl -format vcf4 ${vcf_file} -outfile ${intermediate_file}

# annotate variation
perl ${ANNOVAR}/table_annovar.pl -out ${outpath}/${name_file} \
-build hg19 ${intermediate_file} ${ANNOVAR}/humandb/ \
-remove -protocol wgEncodeGencodeBasicV19,dbnsfp35c -operation gx,f -nastring .

rm ${intermediate_file}

# grep only coding
cat ${outpath}/${name_file}".hg19_multianno.txt" | grep -wv intergenic | grep -wv ncRNA_intronic| grep -vw ncRNA_exonic | grep -wv ncRNA | grep -wv intronic | grep -v UTR | grep -wv downstream | grep -wv upstream | gzip > ${finalfile}
rm ${outpath}/${name_file}".hg19_multianno.txt"