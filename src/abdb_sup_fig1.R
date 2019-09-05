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

vgene.distribution <- function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2_imgt_vgene.csv'
  df = read_csv(infile)
  print(head(df))
  df[df == 'Mus musculus_C57BL/6'] = 'Mus musculus (C57BL/6)'
  countdf = count(df, pdbid, imgt_species, imgt_vgenename)
  print(countdf)
  countdfspec = count(df, imgt_vgenename, imgt_species)
  countdfspec2 = count(countdfspec, imgt_species)
  countdfspec2$x = c(12,20)
  countdfspec2$y = c(60, 40)
  print(countdfspec2)
  total_labels = sprintf('Total number of V genes: %s', countdfspec2$nn)
  ggplot(data=countdf) + 
  geom_bar(mapping = aes(x=imgt_vgenename, fill = imgt_species)) + 
  # coord_polar() + 
  facet_wrap(~ imgt_species, ncol=1, nrow=2, scales='free') +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 90, size = 5)) + 
  scale_fill_manual(values = c(alpha(abcolor3,0.5), alpha(abcolor3))) +
  labs(x= 'IMGT V gene', y='V gene usage (count)') +
  geom_text(data=countdfspec2, mapping = aes(x=x, y=y, label = total_labels)) +
  theme(legend.position = 0)
  outpdf(infile, 'usage_boxplot', width = 10, height=16)

}


absegment_sequence_diversity <- function() {
  infile = 'abdb_outfiles_2019/abdb_segment_absequence_full_vgene_imgt_vgene_ld.csv'
  df = read.csv(infile)
  df = df[seq(100), ]
  df <- df %>% group_by(segment) %>% mutate(median = median(ld))
  print(head(df))
  meds = aggregate(df[,6], list(segment=df$segment), median)
  print(meds)
  med_labels = paste0('Median: ', meds$ld)
  print(med_labels)
  # stop()
  ggplot(data=df) +
  geom_bar(stat='count', mapping = aes(x=ld, fill=segment), show.legend = TRUE) + 
  facet_wrap(~ segment, ncol = 7) + 
  geom_text(data=meds, mapping=aes(x=15, y=10, label = med_labels)) + 
  theme(axis.text.y = element_blank()) + 
  scale_fill_manual(values=my_spectral) +
  coord_polar() + 
  scale_y_log10() +
  theme(legend.position = 0) + 
  # theme(panel.grid.major = element_line(colour = 'grey')) +
  labs(x= ' ', y='Levenshtein distance', fill = 'Region')
  outpdf (infile, 'polar_bar', width = 12, height = 10)
}


ag_sequence_diversity <- function() {
  infile = 'abdb_outfiles_2019/heavy_light_ag_aaseq_antigen_full_ld.csv'
  df = read_csv(infile)
  # df = df[seq(100), ]
  print(df)
  med = median(df$ld)
  label = sprintf('Median: %s', med)
  ggplot(data=df) +
  geom_bar(mapping = aes(x=ld)) + 
  coord_polar() + 
  scale_y_log10() +
  theme(legend.position = 0) + 
  theme(axis.text.y = element_blank()) + 
  # theme(panel.grid.major = element_line(colour = 'grey')) +
  geom_text(mapping = aes(x=1, y=1, label = label), color='Red') + 
  labs(x= ' ', y='Levenshtein distance', fill = 'Region')
  outpng(infile, 'ag_polar_bar', width = 10, height = 10)
}


agsegment_sequence_diversity <- function() {
  infile = 'abdb_outfiles_2019/heavy_light_ag_aaseq_antigen_epitope_ld.csv'
  df = read.csv(infile)
  # df = df[seq(100), ]
  df <- df %>% group_by(region) %>% mutate(median = median(ld))
  print(head(df))
  meds = aggregate(df[,6], list(region=df$region), median)
  print(meds)
  med_labels = paste0('Median: ', meds$ld)
  print(med_labels)
  # stop()
  ggplot(data=df) +
  geom_bar(stat='count', mapping = aes(x=ld, fill=region), show.legend = TRUE) + 
  facet_wrap(~ region, ncol = 2) + 
  geom_text(data=meds, mapping=aes(x=15, y=10, label = med_labels)) + 
  theme(axis.text.y = element_blank()) + 
  scale_fill_manual(values=my_spectral) +
  coord_polar() + 
  scale_y_log10() +
  theme(legend.position = 0) + 
  # theme(panel.grid.major = element_line(colour = 'grey')) +
  labs(x= ' ', y='Levenshtein distance', fill = 'Region')
  outpdf (infile, 'polar_bar', width = 6, height = 12)
}


# run stuff
vgene.distribution()
# absegment_sequence_diversity()
# ag_sequence_diversity()
# agsegment_sequence_diversity()



