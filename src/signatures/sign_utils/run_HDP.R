args = commandArgs(trailingOnly = TRUE)
# we will run it 14 times

i <- as.numeric(args[1])

library(hdp)

samples <- read.table("data/hartwig/signatures/matrix/Colon-Rectum.HDP.txt", sep ='\t',  header = T, stringsAsFactors =F, row.names = 1)

sigs <- read.csv("data/signatures_files/PCAWG/sigProfiler_SBS_signatures_2018_03_28.csv", header = T, row.names = 1, stringsAsFactors = F)

drops <- c("Type","SubType")
sigs <- sigs[ , !(names(sigs) %in% drops)]

# prior with SBS1s
gdsigs <- c("SBS1", "SBS1")

prior_sigs <- as.matrix(sigs[,gdsigs])

nps <- ncol(prior_sigs)

nsamples <- nrow(samples)

ttype_prior <- hdp_prior_init(prior_distn = prior_sigs,
                              prior_pseudoc = rep(1000, nps),
                              hh=rep(1, 96),
                              alphaa=c(1, 1),
                              alphab=c(1, 1))

ttype_prior <- hdp_addconparam(ttype_prior,
                               alphaa = c(1,1),
                               alphab = c(1,1))

ttype_prior <- hdp_adddp(ttype_prior,
                         numdp = nsamples+1,
                         ppindex = c(1, rep(1+nps+1, nsamples)),
                         cpindex = c(3, rep(4, nsamples)))

ttype_prior <- hdp_setdata(ttype_prior,
                           dpindex = (1+nps+1)+1:nsamples,
                           samples[1:nsamples,])

ttype_activated <- dp_activate(ttype_prior,
                               dpindex = (1+nps+1)+0:nsamples,
                               initcc = nps+10,
                               seed = i*1000)

test_chlist <- hdp_posterior(ttype_activated,
                             burnin = 5000,
                             n = 75,
                             space = 200,
                             cpiter = 3,
                             seed = i*1e6,
                             verbosity = 1)

assign(paste0("hdp_", i), test_chlist)
save(list=paste0("hdp_", i), file = paste('data/hartwig/signatures/HDP/ColonRectum.', i, '.hdp.RData', sep =''))