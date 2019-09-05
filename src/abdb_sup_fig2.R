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

paratope_residue_distribution = function() {
  infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments.csv'
  df = read_csv(infile)
  print(df)
  df[df == 'H'] = 'Heavy'
  df[df == 'L'] = 'Light'
  countdf = count(df, abres,segment, abchain) %>% group_by(abchain) %>% mutate(percent = n/sum(n)*100)
  print(countdf)
  csumlabeldf = count(df, abres, abchain) %>% group_by(abchain) %>% mutate(percent = round(n/sum(n)*100, 1))
  print(csumlabeldf)
  sums = csumlabeldf %>% group_by(abchain) %>% summarise(sum = sum(percent))
  print(sums)
  # stop()
  ggplot(data=countdf) + 
    geom_bar(mapping = aes(x=abres, y=percent, fill=segment), stat='identity') + 
    facet_wrap(~abchain) + 
    coord_flip() + 
    scale_fill_manual(values = my_spectral) + 
    labs(y = '', x = 'Residue name', fill = 'Region') + 
    theme(legend.position = 0) +
    geom_text(data=csumlabeldf, mapping = aes(x=abres, y=percent + 0.9, label = percent))
    # scale_y_continuous(labels = percent)
  outpdf(infile, 'paratope_residue_distribution', width = 11)
}


epitope_residue_distribution = function() {
  infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments.csv'
  df = read_csv(infile)
  df[df == 'H'] = 'Heavy'
  df[df == 'L'] = 'Light'
  print(df)
  countdf = count(df, agres,segment, abchain) %>% group_by(abchain) %>% mutate(percent = n/sum(n)*100)
  print(countdf)
  csumlabeldf = count(df, agres, abchain) %>% group_by(abchain) %>% mutate(percent = round(n/sum(n)*100, 1))
  print(csumlabeldf)
  sums = csumlabeldf %>% group_by(abchain) %>% summarise(sum = sum(percent))
  print(sums)
  ggplot(data=countdf) + 
    geom_bar(mapping = aes(x=agres, y=percent, fill=segment), stat='identity') + 
    facet_wrap(~abchain) + 
    coord_flip() +
    scale_fill_manual(values = my_spectral) + 
    labs(y = '', x = 'Residue name', fill='Region') + 
    theme(legend.position = 0) +
    geom_text(data=csumlabeldf, mapping = aes(x=agres, y=percent + 0.9, label = percent))
    # scale_y_continuous(labels = percent)
  outpdf(infile, 'epitope_residue_distribution', width = 11)
}


epitope_gap_residue_distribution = function() {
  infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_len_gap_residue.csv'
  df = read_csv(infile)
  print(df)
  df[df == 'H'] = 'Heavy'
  df[df == 'L'] = 'Light'
  countdf = count(df, residue,segment, abchain) %>% group_by(abchain) %>% mutate(percent = n/sum(n)*100)
  print(countdf)
  csumlabeldf = count(df, residue, abchain) %>% group_by(abchain) %>% mutate(percent = round(n/sum(n)*100, 1))
  print(csumlabeldf)
  sums = csumlabeldf %>% group_by(abchain) %>% summarise(sum = sum(percent))
  print(sums)
  ggplot(data=countdf) + 
    geom_bar(mapping = aes(x=residue, y=percent, fill=segment), stat='identity') + 
    facet_wrap(~abchain) + 
    coord_flip() +
    scale_fill_manual(values = my_spectral) + 
    labs(y = 'Percentage (%)', x = 'Residue name', fill='Region') +
    geom_text(data=csumlabeldf, mapping = aes(x=residue, y=percent + 0.9, label = percent))
    # scale_y_continuous(labels = percent)
  outpdf(infile, 'epitope_gap_residue_distribution', width = 11)
}


paratope_gap_residue_distribution = function() {
  infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_len_gap_residue.csv'
  df = read_csv(infile)
  print(df)
  df[df == 'H'] = 'Heavy'
  df[df == 'L'] = 'Light'
  countdf = count(df, residue,segment, abchain) %>% group_by(abchain) %>% mutate(percent = n/sum(n)*100)
  print(countdf)
  csumlabeldf = count(df, residue, abchain) %>% group_by(abchain) %>% mutate(percent = round(n/sum(n)*100, 1))
  print(csumlabeldf)
  sums = csumlabeldf %>% group_by(abchain) %>% summarise(sum = sum(percent))
  print(sums)
  ggplot(data=countdf) + 
    geom_bar(mapping = aes(x=residue, y=percent, fill=segment), stat='identity') + 
    facet_wrap(~abchain) + 
    coord_flip() +
    scale_fill_manual(values = my_spectral) + 
    labs(y = 'Percentage (%)', x = 'Residue name', fill = 'Region') + 
    geom_text(data=csumlabeldf, mapping = aes(x=residue, y=percent + 0.9, label = percent))
    # scale_y_continuous(labels = percent)
  outpdf(infile, 'epitope_gap_residue_distribution', width = 11)
}


#run stuff
paratope_residue_distribution()
epitope_residue_distribution()
epitope_gap_residue_distribution()
paratope_gap_residue_distribution()


