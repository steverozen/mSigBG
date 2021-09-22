library(tidyr)
library(lsa)
library(ggsignif)

data.dir <- file.path('data-raw')
source(paste0(data.dir,'/syn.data.function.R'))

#target: 'SBSXX' form
#bg: 'mSigBG::background.info$HepG2', spectra form
#ratio: num of target muts/num of muts in bg, must be <= 1, default = 1
#num: number of synthesized cell lines
target='SBS22'
bg = mSigBG::background.info$HepG2
name.of.bg = strsplit(colnames(bg$background.sig), '[.]')[[1]][1]
ratio = 1
num = 10

list1 <- generate.syn.spec(num, ratio = ratio, target, bg)
#output: list(output.syn.spec, diff.muts.in.bg, diff.muts.in.target)

ret<- 
  SeparateSignatureFromBackground(
    spectra     = list1[[1]][,-1],
    bg.sig.info = bg,
    m.opts      = NULL,
    start.b.fraction = 0.5)

inferred.sig <- 
  ICAMS::as.catalog(
    object         = ret$inferred.target.sig, 
    catalog.type   = "counts.signature",
    region         = "genome",
    abundance      = attr(bg$background.sig,
                          "abundance"),
    infer.rownames = TRUE)

#draw inferred target signature spectrum
colnames(inferred.sig) <- paste0("Inferred-", target, "-sig")
two.sigs <- cbind(inferred.sig, bg$background.sig)
ICAMS::PlotCatalogToPdf(two.sigs,file=paste0('inferred.',ratio, target, '.sig.in.',name.of.bg, '.pdf'))

#draw inferred bg and target counted spectrum
inferred.bg.spectra <- 
  round(bg$background.sig %*% matrix(ret$exposures.to.bg.sig, nrow = 1))
inferred.bg.spectra <- ICAMS::as.catalog(inferred.bg.spectra,
                                         catalog.type   = "counts",
                                         region         = "genome",
                                         abundance = attr(ideal.syn.spec, "abundance"),
                                         infer.rownames = TRUE)
colnames(inferred.bg.spectra) <- paste0(name.of.bg, '.', sub('.*([1-9][0-9]?).*','cl\\1',colnames(list1[[1]][,-1])), "-inf-bg-spect")

inferred.target.spectra <- ideal.syn.spec[,rep(1,time = num)] - inferred.bg.spectra
colnames(inferred.target.spectra) <- paste0(name.of.bg,'.', sub('.*([1-9][0-9]?).*','cl\\1',colnames(list1[[1]][,-1])), "-inf-target-spect")
ICAMS::PlotCatalogToPdf(cbind(inferred.bg.spectra,inferred.target.spectra),file = paste0(name.of.bg, '_inferred_count_with_', ratio, target,'.pdf'))

#draw scatter plots for count diff between inferred num of muts and ideal num of muts
inferred.muts <- cbind(ret$exposures.to.target.sig, ret$exposures.to.bg.sig, list1[[3]], list1[[2]])
colnames(inferred.muts) <- c(paste0('infer.target.', target), paste0('infer.bg.', name.of.bg),
                             paste0('assigned.target.', target), paste0('assigned.bg.', name.of.bg))
inferred.muts <- data.frame(names = sub('.*([1-9][0-9]?).*','cl\\1',rownames(inferred.muts)),inferred.muts)
rownames(inferred.muts) <- NULL;inferred.muts
inferred.muts$sum.inferred <- rowSums(cbind(inferred.muts$infer.target, inferred.muts$infer.bg))
inferred.muts$sum.ideal <- rowSums(cbind(inferred.muts$assigned.target, inferred.muts$assigned.bg))

forscatterplot = inferred.muts[,1:5]%>% #dataframe row to column transformation
  pivot_longer(cols = !1,
               names_to = "source",
               values_to = "count")
title = paste0('Diff between inferred count and ideal count for ', ratio, target, ' in ', name.of.bg)
forscatterplot$names <- factor(forscatterplot$names , levels = unique(forscatterplot$names)) #stop ggplot2 from reordering x axis
my_comparisons <- list(c("assigned.bg.HepG2", "infer.bg.HepG2"), c("assigned.target.SBS22", "infer.target.SBS22"))
p <- ggplot(data = forscatterplot) +
  geom_point(mapping = aes(x=names,y=count,color=source),stat="identity")+
  scale_color_manual(values = c(infer.target.SBS22 = '#3399FF', infer.bg.HepG2= '#FF9999', 
                                assigned.target.SBS22 = '#009999', assigned.bg.HepG2 = '#CC0000'))+
  #stat_compare_means(comparisons = my_comparisons)+
  #stat_compare_means(label.y = 50)+
  ggtitle(title)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = p, filename = paste0(title, '.jpg'))
p

forscatterplot$source <- factor(forscatterplot$source, levels=
                                  c('assigned.bg.HepG2', 'infer.bg.HepG2',
                                    'assigned.target.SBS22', 'infer.target.SBS22'))
p2 <- ggplot(data = forscatterplot,mapping = aes(x=source, y=count,fill = source)) +
  geom_boxplot()+
  geom_point()+
  geom_signif(comparisons = my_comparisons)+
  ggtitle(title)+
  theme_bw()+
  scale_fill_manual(values = c(infer.target.SBS22 = '#3399FF', infer.bg.HepG2= '#FF9999', 
                               assigned.target.SBS22 = '#009999', assigned.bg.HepG2 = '#CC0000'))
ggsave(plot = p2, filename = paste0(title, 'box.jpg'))
p2
#cosine distance between target signatures
cosine.dis.target <- lsa::cosine(inferred.sig[,1],target.sig[,1])
cosine.dis.bg <- lsa::cosine(bg$background.sig[,1], 
                             mSigBG::MeanOfSpectraAsSig(inferred.bg.spectra)$mean.sig[,1])

print(as.data.frame(x =c(cosine.dis.bg,cosine.dis.target),row.names=c('bg', 'target'),col.names = 'cosine distance'))
