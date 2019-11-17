# Copyright (C) 2019 Institute for Research in Biomedicine,  Moritz Gerstung
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


library(VariantAnnotation)
library('boot')
library('RColorBrewer')
args = commandArgs(trailingOnly=TRUE)


prepare_vcf <- function(c_muts){
  suppressWarnings(
    vcf <- VCF(rowRanges = GRanges(
      seqnames = c_muts$CHROM,
      ranges = IRanges(
        start = c_muts$POS,
        width = as.integer(rep(1, nrow(c_muts)))
      ),
      strand = rep("+", nrow(c_muts)),
      REF = c_muts$REF,
      ALT = c_muts$ALT),
      info = DataFrame(
        t_ref_count = c_muts$REF_COUNTS,
        t_alt_count = c_muts$VAR_COUNTS,
        gene = rep("NA", nrow(c_muts)),
        effect = rep("NA", nrow(c_muts)),
        mut = "SNV",
        state = rep("NA", nrow(c_muts))

      )
    )
  )
}

#---------------
# EXTRA FILES
#--------------

# run the script
load('src/timing/refLengths_chrOffset.RData')
chrOffset <- read.csv("src/timing/chromOff.tsv", sep ='\t', header = TRUE, check.names=FALSE)

names <- colnames(chrOffset)
chrOffset <- as.numeric(chrOffset)
names(chrOffset)<-names


#---------------
# LOADING FILES
#---------------

# MUTATION FILE
mutation_file <-args[1]

# name of the sample
name_sample <- as.character(strsplit(basename(mutation_file), split = '.', fixed = TRUE)[[1]][1])

# CNV FILE, output of cnas_ccf script
c_cnv_file = paste('data/hartwig/clonal_structure/', name_sample, '.cna.gz', sep ='')

# CLONAL RECONSTRUCTION FILE, output of clustering_ccf script
clonal_reconstruction_file = paste('data/hartwig/subclonal/', name_sample, '.subclonal.txt', sep ='')

# check if there was a successful clonal reconstruction (not in blacklisted folder)
if (file.exists(clonal_reconstruction_file)) {

  # create outpath if it doesn't exist
  dir.create(file.path('data/hartwig/timing/'), showWarnings = FALSE)

  #------------
  # LOAD CNVs
  #------------
  print("Loading CNVs...")
  bb <- read.table(c_cnv_file, sep ='\t', header = TRUE)
  bb$chromosome<- gsub( "chr", "", as.character(bb$chromosome))  # get the gender

  # get the gender info
  gender <- as.character(bb$gender[1])

  bb <- GRanges(
    seqnames = bb$chromosome,
    ranges = IRanges(
      start = bb$start,
      end = bb$end
    ),
    mcols = bb[,c("major_cn", "minor_cn", "clonal_frequency", "n.snv_mnv")]
  )
  colnames(mcols(bb)) <- c("major_cn", "minor_cn", "clonal_frequency", "n.snv_mnv")
  bb$total_cn <- bb$major_cn + bb$minor_cn

  #---------------
  # LOAD CLUSTERS
  #---------------
  print("Loading clusters...")
  c_clusters <- read.table(clonal_reconstruction_file, sep ='\t',
                           header = TRUE)
  #----------------
  # LOAD MUTATIONS
  #----------------
  print("Loading mutations...")
  c_muts <- read.table(mutation_file, sep ='\t', header = TRUE)
  c_muts$CHROM<- gsub( "chr", "", as.character(c_muts$CHROM))# get the gender
  vcf<-prepare_vcf(c_muts)


  #---------------
  # Get WGD STATUS
  #---------------
  # Hom, Ploidy, WGstatus
  hom <- sum(width(bb) * (bb$minor_cn == 0) * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
  ploidy <- sum(width(bb) * (bb$major_cn + bb$minor_cn) * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
  isWgd <- 2.9 - 2*hom <= ploidy

  # load MutationTimeR slightly modified functions
  source('src/timing/MutationTime.R')

  #----------------------
  # MUTATIONAL ANALYSIS
  #----------------------
  print("Compute Mutational Time....")

  # compute mutational timing
  MT <- computeMutCn(vcf, bb,c_clusters, c_clusters[1,1], isWgd = isWgd, gender = gender)

  print("Classify Mutations...")
  MT$D$CLS <- classifyMutations(as.data.frame(MT$D))

  c_muts$timing_class <- MT$D$CLS
  c_muts$pMutCN <- MT$D$pMutCN
  c_muts$pGain <- MT$D$pGain
  c_muts$pSingle <- MT$D$pSub
  c_muts$pSub <- MT$D$pSub
  c_muts$pMutCNTail <- MT$D$pMutCNTail
  c_muts$WGD <- isWgd
  c_muts$ploidy <- ploidy
  c_muts$hom <- hom

  name_outfile = (paste('data/hartwig/timing/', name_sample, '.mutationaltiming.tsv.gz', sep =''))
  write.table(c_muts, gzfile(name_outfile), sep ='\t', col.names =  TRUE, row.names = FALSE, quote = FALSE)

  #-----------------
  # QQPLOT ANALYSIS
  #-----------------

  png(filename = (paste('data/hartwig/timing/', name_sample, '.qqplot.png', sep ='')), width=300, height=300)
  par(mfrow = c(1,1))
  qqnorm(qnorm(MT$D$pMutCNTail), xlim=c(-5,5), ylim=c(-5,5), pch=16)
  abline(0,1, col='red')
  dev.off()

  pval_array <- MT$D$pMutCNTail

  reduced <- pval_array[(pval_array>0.15)&(pval_array<0.85)]
  uniform <- runif(10000, min = 0.15, max = 0.85)

  pl <- ks.test(reduced, uniform )

  qqplot_out <- data.frame("D"=pl$statistic, "p-value"=pl$p.value)
  name_outfile = (paste('data/hartwig/timing/', name_sample, '.qqplot.tsv', sep =''))
  write.table(qqplot_out, name_outfile, sep ='\t', col.names =  TRUE, row.names = FALSE, quote = FALSE)

  #----------------------
  # PLOTTING WGD TIMING
  #----------------------

  n.snv_mnv <- countOverlaps(bb,vcf)
  info(vcf) <- cbind(info(vcf), MT$D)
  info(vcf)$CLS<-MT$D$CLS
  time <- bbToTime(bb, timing_param=MT$P, n.snv_mnv=n.snv_mnv)

  bb$n.snv_mnv <- n.snv_mnv
  bb$time <- time$time
  bb$time.lo <- time$time.lo
  bb$time.up <- time$time.up
  bb$time.2nd <- time$time.2nd
  bb$time.2nd.lo <-time$time.2nd.lo
  bb$time.2nd.up <-time$time.2nd.up
  bb$time.star <-time$time.star
  bb$time.type <-time$type

  name_outfile = (paste('data/hartwig/timing/', name_sample, '.mutationaltiming.bb.tsv.gz', sep =''))
  write.table(bb, gzfile(name_outfile), sep ='\t', col.names =  TRUE, row.names = FALSE, quote = FALSE)

  png(filename = (paste('data/hartwig/timing/', name_sample, '.cnaplot.png', sep ='')),
      width=2500, height=2500, res = 300)

  par(mfrow=c(1,1), xpd=TRUE)# par()
  molecular_timing <-plotSample(vcf, bb, ploidy=ploidy, hom = hom, ID = name_sample, IS_WGD = isWgd, clusters = c_clusters)
  dev.off()

  #==========================
  # SAVING MOLECULAR TIMING
  #==========================

  mol_out <- data.frame("molecular_timing"=molecular_timing, "sample"=name_sample, 'isWGD'=isWgd)
  name_outfile = (paste('data/hartwig/timing/', name_sample, '.WGD_time.tsv', sep =''))
  write.table(mol_out, name_outfile, sep ='\t', col.names =  TRUE, row.names = FALSE, quote = FALSE)
}
