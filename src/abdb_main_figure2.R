setwd('/Users/rahmadakbar/greifflab/aims/aimugen/src')
# detach('package:themeakbar', unload = TRUE)
library('themeakbar')
library('janitor')
library('ggplot2')
library('pheatmap')
library('beeswarm')
library(data.table)
library('dplyr')
library('rlist')
library('scales')
library('reshape2')
library('VennDiagram')
library('tidyr')
library('RCy3')
library(RColorBrewer)
library('igraph')
library('tidyverse')
library('viridis')
library('poweRlaw')
library('fossil')
library('gtools')
library('vegan')
library('iNEXT')



theme_set(themeakbar())
my_spectral <- colorRampPalette(brewer.pal(8,'Spectral'))(14)
abcolor = '#6464FF'
abcolor2 = '#1010FF'
abcolor3 = '#3399CC'
agcolor = '#C0C0C0'
blues = colorRampPalette(brewer.pal(8,'Blues'))(14)[2:10]
oranges = brewer.pal(6,'Oranges')
greys = brewer.pal(6,'Greys')
blues2 = brewer.pal(6,'Blues')
blues3 = c('#9010FF', '#1010FF', '#C8C8FF')
blues4 = c(blues3, blues2[2:6])
cdrfr_colors = c(oranges,blues)
rewcolors = c(greys[3:5], c('#9010FF', '#1010FF', '#C8C8FF'))
rewcolors2  = c(greys[3:4], blues3[1:2])
blues5 = c('#1010FF','#1010FF','#1010FF',greys[3],'#1010FF')
outfigdir = 'abdb_figures_2019'


motif.length.distribution <- function(){
  infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_len.csv'
  df = read.csv(infile)
  print(head(df))
  print(group_by(df, segment))
  # df$motif_len = as.numeric(df$motif_len)
  # df <- df %>% group_by(segment) %>% mutate(median = median(motif_len))
  # meds = aggregate(df[,6], list(segment=df$segment), median)
  # print(meds)
  # med_labels = paste0('median: ', meds$ld)
  # print(med_labels)
  ggplot(data=df) + 
    # geom_bar(mapping = aes(x=motif_len), fill=my_spectral[1]) + 
    # geom_bar(stat='count', mapping = aes(x=motif_len), fill='red') + 
    geom_density(mapping = aes(x=motif_len), fill=alpha(agcolor,0.7), color=alpha(agcolor, 0.7)) + 
    facet_wrap(~ segment, ncol = 3, scales = 'free') + 
    labs(x='Motif length (# of residues + # of gaps)', y='Density') + 
    coord_flip()
  outname = 'abdb_figures_2019/motif_length_distribution_epitope.pdf'
  ggsave(outname, height = 9, width=7)
  system(sprintf('open %s', outname))
}
  

ab_ag_motiflen_correlation  = function(){
  # examine ab ag angle per region
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged.csv'
  df = read.csv(infile)
  print(head(df))
  df = drop_na(df)
  segments = sort(unique(df$segment))
  print(segments[1])
  cors_data = sapply(1:length(segments), function(x) cor.test(subset(df, df$segment == segments[x])$ab_motiflen,subset(df, df$segment == segments[x])$ag_motiflen))
  cors = sapply(1:dim(cors_data)[2], function(x) round(cors_data[,x]$estimate,2))
  print(cors_data[,1])
  pvals = sapply(1:dim(cors_data)[2], function(x) signif(cors_data[,x]$p.value, digits = 1))
  # cors_pearson = df %>% group_by(segment) %>% summarize(cor = cor(p_angle, e_angle, method='pearson', use = 'complete.obs'))
  cors_pearson = df %>% group_by(segment) %>% summarize(cor = cor(ab_motiflen, ag_motiflen, method='pearson', use = 'complete.obs'))
  cors_spearman = df %>% group_by(segment) %>% summarize(cor = cor(ab_motiflen, ag_motiflen, method='spearman', use = 'complete.obs'))
  print(cors_pearson)
  print(cors_spearman)
  label_pearson = paste0('Pearson Cor. ', round(cors_pearson$cor, 2))
  label_spearman = paste0('Spearman Cor. ', round(cors_spearman$cor, 2))
  print(label_pearson)
  ggplot(data=df, mapping = aes(x=ab_motiflen, y=ag_motiflen)) + 
    geom_point(mapping = aes(x = ab_motiflen, y= ag_motiflen), color = alpha(my_spectral[3],0.5), size=0.5) + 
    geom_text(data = cors_pearson, mapping = aes(x=10, y=20, label=label_pearson)) +
    geom_text(data = cors_pearson, mapping = aes(x=10.5, y=17, label=label_spearman)) +
    facet_wrap(~ segment, ncol = 3) + 
    geom_smooth(method = 'lm', mapping = aes(x=ab_motiflen, y= ag_motiflen), color=alpha(abcolor, 0.3)) + 
    geom_density2d(color = alpha(my_spectral[3],0.7)) +
    labs(x='Motif length (paratope)', y='Motif length (epitope)')
  outname = 'abdb_figures/ab_ag_motiflen_correlation.pdf'
  ggsave(outname, height = 9, width=7)
  system(sprintf('open %s', outname))
  
}

ab.ag.venn = function(){
  infile_ab = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx.csv' 
  infile_ag = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx.csv' 
  # infile_ag = 'abdb_outfiles/respairs_epitope_segment_notationx_top-1_gap_patterns.csv' 
  dfab = read.csv(infile_ab, stringsAsFactors = FALSE)
  dfag = read.csv(infile_ag, stringsAsFactors = FALSE)
  motif_ab = dfab$gap_patterns
  print(motif_ab)
  motif_ag = dfag$gap_patterns
  x = list(motif_ab,motif_ag)
  #random sampling
  motif_ab_random = sample(motif_ab, 391)
  print(motif_ab_random)
  sapply(1:10, function(xet) get.venn.partitions(x))
  ###
  # outname ='abdb_figures/ab_ag_venn_random.tiff'
  inname = strsplit(tail(strsplit(infile_ab, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_venn.tiff', outfigdir,inname)
  # print(outname)
  names(x) = c('motif_ab', 'motif_ag')
  venn.plot <- venn.diagram(x=x, filename=outname,
                            output = FALSE ,
                            imagetype="tiff" ,
                            # height = 800 , 
                            # width = 900 , 
                            resolution = 300,
                            compression = "lzw",
                            lwd = 2,
                            lty = 'blank',
                            fill = c(abcolor3, agcolor),
                            alpha=c(0.7,0.7),
                            rotation.degree = 180,
                            cex = 5, 
                            cat.fontface=3,
                            cat.cex = 3,
                            # cat.dist = c(-0.05,-0.15),
                            cat.just = list(c(1.4, 0), c(-1.4, -3.4)),
                            category.names=c("Paratope", "Epitope")
                            )
  system(sprintf('open %s', outname))
  
}

motif_coverage = function(){
  infile = 'abdb_outfiles_2019/motif_coverage_skip30.csv' 
  df = read.csv(infile, stringsAsFactors = FALSE)
  df = transform(df, unique_motif = as.factor(unique_motif))
  print(sapply(df, class))
  df$source[df$source == 'ab'] <- 'Paratope'
  df$source[df$source == 'ag'] = 'Epitope'
  df$source = factor(df$source, levels = c("Paratope", 'Epitope')) 
  print(head(df))
  medpara = aggregate(df[df$source=='Paratope',4], by= list(df[df$source=='Paratope',]$pdbidx),median)
  medepi = aggregate(df[df$source=='Epitope',4], by= list(df[df$source=='Epitope',]$pdbidx),median)
  print(tail(medpara))
  print(tail(medepi))
  ggplot(df, aes(x=pdbidx, y=coverage)) + 
    # geom_point()+
    # geom_point(mapping = aes(x=pdbidx, y=coverage)) + 
    geom_boxplot(mapping = aes(x= reorder(strpdbidx, pdbidx), y=coverage*100, fill= source), lwd=0, outlier.shape = NA, flatten=0) +
    facet_wrap(~source, scales = 'free') +
    scale_fill_manual(values = c(abcolor3, agcolor), name = "Motif source", labels = c('Paratope', 'Epitope')) + 
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = 'Number of structures', y = 'Interaction motif coverage')
    # coord_flip()
    # geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE)
  outname = 'abdb_figures_2019/motif_coverage.pdf'
  ggsave(outname, height = 10, width = 15)
  system(sprintf('open %s', outname))
}

chao_diversity = function(){
  #compute diversity estimate using chao2
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged.csv'
  df = read_csv(infile)
  print(df)
  dfab = df %>% group_by(ab_motif) %>%tally()
  dfab = dfab[order(-dfab$n),]
  chao_ab = chao1(dfab$n)
  print(dfab)
  dfag = df %>% group_by(ag_motif) %>%tally()
  dfag = dfag[order(-dfag$n),]
  i.out = iNEXT(dfab$n, datatype = 'abundance')
  chao_abinext = ChaoRichness(dfab$n)
  print(chao_ab)
  print(chao_abinext)
  ggiNEXT(i.out,type=1)
  chao_ag = chao1(dfag$n)
  print(chao_ag)
  print(dfag)
  obs_ab = nrow(dfab)
  obs_ag = nrow(dfag)
  divdf = data.frame('divs' = c(obs_ab, obs_ag, chao_ab, chao_ag), 'type' = c('Observed', 'Observed', 'Estimated', 'Estimated'), 'source' = 
                       factor(c('Paratope', 'Epitope', 'Paratope', 'Epitope'), levels = c('Paratope', 'Epitope')))
  print(divdf)
  ggplot(data = divdf) +
    geom_bar(stat = 'identity', mapping = aes(reorder(type,divs), divs, fill = source), position = 'dodge') + 
    facet_wrap(~  source, scale = 'free') + 
    scale_fill_manual(values = c(abcolor3, agcolor)) + 
    geom_text(mapping = aes(reorder(type, divs), y=c(410,1750,750, 10650) ), label = round(divdf$divs)) + 
    labs(x = '', y = 'Interaction motif diveristy', fill = 'Motif source')
  outname = 'abdb_figures_2019/chao_diversity.pdf'
  ggsave(outname, height = 10, width = 8)
  system(sprintf('open %s', outname))
  
}



# runstuff
# motif.length.distribution()
# ab_ag_motiflen_correlation()
# ab.ag.venn() 
# motif_coverage()
chao_diversity()
  
