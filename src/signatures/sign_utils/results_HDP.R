library(hdp)
# this is partially based on Nicola Robert's examples.


chlist <- vector("list", 14)
for (i in 1:14 ) {
  print(i)

  load(paste('data/hartwig/signatures/HDP/ColonRectum.', i, '.hdp.RData', sep =''))
}
chlist[[1]] <- hdp_1
chlist[[2]]<- hdp_2
chlist[[3]]<- hdp_3
chlist[[4]]<- hdp_4
chlist[[5]]<- hdp_5
chlist[[6]]<- hdp_6
chlist[[7]]<- hdp_7
chlist[[8]]<- hdp_8
chlist[[9]]<- hdp_9
chlist[[10]]<- hdp_10
chlist[[11]]<- hdp_11
chlist[[12]]<- hdp_12
chlist[[13]]<- hdp_13
chlist[[14]]<- hdp_14

mut_example_multi <- hdp_multi_chain(chlist)
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")

hdpsample <- hdp_extract_components(mut_example_multi)

par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(hdpsample, bty="L")
trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))
# pick your colours
mut_colours <- c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70')

plot_comp_distn(hdpsample, cat_names=trinuc_context,
                grouping=group_factor, col=mut_colours,
                col_nonsig="grey80", show_group_labels=TRUE)

# modified to include number of samples
dpindices=4+(1:488)

plot_dp_comp_exposure(hdpsample, main_text="Colon-Rectum",
                      dpindices=dpindices,
                      col= c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(12, "Set3")),
                      incl_nonsig=FALSE,
                      ylab_numdata = 'SNV count', ylab_exp = 'Signature exposure',
                      leg.title = 'Signature')

# save exposures
main_text = NULL
dpnames = NULL

incl_nonsig = FALSE
incl_comp0 = TRUE
dp_distn <- comp_dp_distn(hdpsample)
ndp <- nrow(dp_distn$mean)
ncomp <- ncol(dp_distn$mean)
incl_numdata_plot = TRUE
if (!is.numeric(dpindices) | any(dpindices%%1 != 0) | any(dpindices <
                                                          1) | any(dpindices > ndp)) {
  stop(paste("dpindices must be integers between 1 and",
             ndp))
}
if (!class(dpnames) %in% c("character", "NULL") | !length(dpnames) %in%
    c(length(dpindices), 0)) {
  stop("dpnames must be a character vector with\n         same length as dpindices, or NULL")
}
if (!class(main_text) %in% c("character", "NULL") | !length(main_text) %in%
    c(1, 0)) {
  stop("main_text must be a character string, or NULL")
}
if (class(incl_numdata_plot) != "logical") {
  stop("incl_numdata_plot must be TRUE or FALSE")
}
if (class(incl_nonsig) != "logical")
  stop("incl_nonsig must be TRUE or FALSE")
par_old <- par(no.readonly = TRUE)
on.exit(par(par_old), add = TRUE)

if (class(hdpsample) == "hdpSampleChain") {
  dps <- dp(final_hdpState(hdpsample))[dpindices]
  pps <- ppindex(final_hdpState(hdpsample))[dpindices]
}
if (class(hdpsample) == "hdpSampleMulti") {
  dps <- dp(final_hdpState(chains(hdpsample)[[1]]))[dpindices]
  pps <- ppindex(final_hdpState(chains(hdpsample)[[1]]))[dpindices]
}

numdata <- sapply(dps, function(x) x@numdata)
dp_order <- order(numdata, decreasing = TRUE)
if (incl_numdata_plot & any(numdata == 0)) {
  stop("Can't have incl_numdata_plot TRUE if\n         one or more dpindices have no associated data item/s")
}
if (length(unique(pps)) > 1) {
  warning("some dpindices have different parent nodes,\n            separate plots may be better")
}
exposures <- t(dp_distn$mean[dpindices, , drop = FALSE])
if (!incl_nonsig) {
  cis <- dp_distn$cred.int[dpindices]
  nonsig <- lapply(cis, function(x) which(x[1, ] == 0))
  for (i in 1:length(nonsig)) {
    exposures[nonsig[[i]], i] <- 0
  }
}
if (!incl_comp0) {
  exposures <- exposures[-1, ]
}

write.table(exposures, file = 'data/hartwig/signatures/HDP/exp_method_HDP.tsv', sep ="\t")

## get processes ##
comp_distn <- comp_categ_distn(hdpsample)
ncomp <- nrow(comp_distn$mean)-1
comp_to_plot <- rownames(comp_distn$mean)

processes_df <- data.frame()
for (ii in seq_along(comp_to_plot)){

  cname <- comp_to_plot[ii]

  # mean categorical distribution (sig), and credibility interval
  sig <- comp_distn$mean[cname,]
  processes_df <- rbind(processes_df,sig)
}
colnames(processes_df) <- 1:96
write.table(processes_df, file = 'data/hartwig/signatures/HDP/processes_method_HDP.tsv', sep ="\t")

