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
# library('fossil')
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

outpdf = function(infile, tag, width=8, height=8){
  # assumes infile has extensins (.csv etc)
  # uses the first bit of the file as name and add the tag
  # uses ggsave, opens pdf plot after saving
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_%s.pdf', outfigdir,inname,tag)
  ggsave(outname, width=width, height = height)
  system(sprintf('open %s', outname))
  print(sprintf('opening %s', outname))
  
}

outpng = function(infile, tag, width=8, height=8){
  # assumes infile has extensins (.csv etc)
  # uses the first bit of the file as name and add the tag
  # uses ggsave, opens pdf plot after saving
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_%s.png', outfigdir,inname,tag)
  ggsave(outname, width=width, height = height)
  system(sprintf('open %s', outname))
  print(sprintf('opening %s', outname))
}



ab_ag_angle = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle.csv'
  df = read_csv(infile)
  df = drop_na(data = df, p_angle, e_angle)
  meds_data = aggregate(df[, c('p_angle', 'e_angle')], list(df$segment), median)
  max_data = aggregate(df[, c('p_angle', 'e_angle')], list(df$segment), median)
  colnames(meds_data)[1] = 'segment'
  df2 = gather(df,'angle_type', 'angle', p_angle, e_angle)
  ggplot(data=df2) + 
  # geom_point(mapping = aes(x=p_angle, y=e_angle))
  geom_density(mapping = aes(x=angle, fill=angle_type, color=angle_type), alpha = 0.7) +
  scale_fill_manual(values= c(agcolor, abcolor3), labels=c('Paratope', 'Epitope')) + 
  scale_color_manual(values= c(agcolor, abcolor3), labels=c('Paratope', 'Epitope')) + 
  facet_wrap(~ segment, ncol = 3, scales = 'free') +
  geom_text(data=meds_data, mapping = aes(x= e_angle, y=0.01, label=paste0('Median: ',e_angle)), size=2, color=alpha('Black', 0.7)) +
  geom_text(data=meds_data, mapping = aes(x= p_angle, y=0.015, label=paste0('Median: ',p_angle)), size=2, color=alpha('Blue', 0.95)) +
  # geom_vline(data=agmeddf, mapping = aes(xintercept=agmeds_data), color=alpha('grey',0.7), linetype='dashed') +
  # geom_vline(data=abmeddf, mapping = aes(xintercept=abmeds_data), color=alpha('blue',0.5), linetype='dashed') +
  # theme(legend.text = element_text(size=14)) + 
  labs(x='Motif angle (degree)', y='Density', fill='Motif source', color='Motif source')
  # outname = 'abdb_figures/angle_distribution.pdf'
  # ggsave(outname, width = 8, height = 8)
  # system(sprintf('open %s', outname))
  outpdf(infile, 'angle_distribution')
}


ab_ag_angle_correlation  = function(){
  # examine ab ag angle per region
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle.csv'
  df = read.csv(infile)
  print(head(df))
  df = drop_na(df)
  segments = sort(unique(df$segment))
  print(segments[1])
  scor_data = df %>% group_by(segment) %>% summarise('correlation'=cor(p_angle, e_angle,method='spearman'))
  pcor_data = df %>% group_by(segment) %>% summarise('correlation'=cor(p_angle, e_angle,method='pearson'))
  pcor_labels = sprintf('Pearson cor. %s', round(pcor_data$correlation,2))
  scor_labels = sprintf('Spearman cor. %s', round(scor_data$correlation,2))
  print(pcor_labels)
  # stop()
  # cors_data = sapply(1:length(segments), function(x) cor.test(subset(df, df$segment == segments[x])$p_angle,subset(df, df$segment == segments[x])$e_angle))
  # cors = sapply(1:dim(cors_data)[2], function(x) round(cors_data[,x]$estimate,2))
  # print(cors_data[,1])
  # pvals = sapply(1:dim(cors_data)[2], function(x) signif(cors_data[,x]$p.value, digits = 1))
  # cors_pearson = df %>% group_by(segment) %>% summarize(cor = cor(p_angle, e_angle, method='pearson', use = 'complete.obs'))
  # print(cors_pearson)
  # stop()
  # label_pearson = paste0('r: ', cors, '; pval: ', pvals)
  # print(label_pearson)
  ggplot(data=df) + 
  geom_point(mapping = aes(x = p_angle, y= e_angle), color = alpha(my_spectral[3],0.5), size=1) + 
  geom_density2d(mapping = aes(x=p_angle, y=e_angle), color = alpha(my_spectral[3])) +
  geom_text(data = pcor_data, mapping = aes(x=10, y=245, label=pcor_labels), hjust=0) +
  geom_text(data = pcor_data, mapping = aes(x=10, y=215, label=scor_labels), hjust=0) +
  facet_wrap(~ segment, ncol = 3, scales = 'free') +
  ylim(0,250) + xlim(1,200) +
  geom_smooth(method = 'lm', mapping = aes(x=p_angle, y= e_angle), color=alpha(abcolor, 0.3)) + 
  labs(x='Paratope angle', y='Epitope angle')
  # outname = 'abdb_figures/angle_correlation.pdf'
  # ggsave(outname)
  # system(sprintf('open %s', outname))
  outpdf(infile, 'angle_correlation')
}

# run stuff
# ab_ag_angle()
ab_ag_angle_correlation()