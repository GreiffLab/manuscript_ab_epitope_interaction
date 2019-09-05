# import stuff
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


species_distribution_antibody = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
  df = read_csv(infile)
  print(df)
  # clean up some messed up species name
  df[df=='LLAMA GLAMA'] = "LAMA GLAMA"
  df[df=='MUS MUSCULUSHOMO SAPIENS'] = "MUS MUSCULUS HOMO SAPIENS"
  df = drop_na(df)
  countdf = count(df, hspecies)
  print(countdf)
  total_label = sprintf('Total species: %s', dim(countdf)[1])
  ggplot(data = countdf) + 
    geom_bar(stat = 'identity', mapping = aes(x=reorder(hspecies, n), y=n), fill=abcolor3) + 
    # theme(axis.text.x = element_text(angle= 90)) + 
    geom_text(mapping = aes(x=hspecies, y=n+90, label = n)) + 
    coord_flip() + 
    labs(x='Species', y='# of paratopes') + 
    geom_text(mapping = aes(x=1, y = 1000, label = total_label))
  
  outpdf(infile, 'hspecies_distribution', width = 10, height = 15)
  
}


species_distribution_antibody_complex = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
  df = read_csv(infile)
  print(df)
  # clean up some messed up species name
  df[df=='LLAMA GLAMA'] = "LAMA GLAMA"
  df[df=='MUS MUSCULUSHOMO SAPIENS'] = "MUS MUSCULUS HOMO SAPIENS"
  df = drop_na(df)
  countdf = count(df, hspecies, pdbid)
  print(countdf)
  countdf2 = count(countdf, hspecies)
  print(countdf2)
  total_label = sprintf('Total species: %s', dim(countdf2)[1])
  ggplot(data = countdf2) + 
    geom_bar(stat = 'identity', mapping = aes(x=reorder(hspecies, nn), y=nn), fill=abcolor3) + 
    # theme(axis.text.x = element_text(angle= 90)) + 
    geom_text(mapping = aes(x=hspecies, y=nn+15, label = nn)) + 
    coord_flip() + 
    labs(x='Species', y='# of structures') + 
    geom_text(mapping = aes(x=1, y = 200, label = total_label))
  
  outpdf(infile, 'hspecies_distribution_complex', width = 10, height = 15)
  
}

species_distribution_antigen = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
  df = read_csv(infile)
  print(df)
  # clean up some messed up species name
  # df[df=='LLAMA GLAMA'] = "LAMA GLAMA"
  # df[df=='MUS MUSCULUSHOMO SAPIENS'] = "MUS MUSCULUS HOMO SAPIENS"
  df = drop_na(df)
  countdf = count(df, ag_species2)
  total_label = sprintf('Total species: %s', dim(countdf)[1])
  print(countdf)
  ggplot(data = countdf) + 
    geom_bar(stat = 'identity', mapping = aes(x=reorder(ag_species2, n), y=n), fill=agcolor) + 
    # theme(axis.text.x = element_text(angle= 90)) + 
    geom_text(mapping = aes(x=ag_species2, y=n+90, label = n)) + 
    coord_flip() + 
    labs(x='  ', y='# of epitopes') + 
    geom_text(mapping = aes(x=3, y = 1000, label = total_label))
  
  outpdf(infile, 'agspecies_distribution', width = 10, height = 15)
  
}


species_distribution_antigen_complex = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
  df = read_csv(infile)
  print(df)
  # clean up some messed up species name
  # df[df=='LLAMA GLAMA'] = "LAMA GLAMA"
  # df[df=='MUS MUSCULUSHOMO SAPIENS'] = "MUS MUSCULUS HOMO SAPIENS"
  df = drop_na(df)
  countdf = count(df, ag_species2, pdbid)
  print(countdf)
  countdf2 = count(countdf, ag_species2)
  total_label = sprintf('Total species: %s', dim(countdf2)[1])
  ggplot(data = countdf2) + 
    geom_bar(stat = 'identity', mapping = aes(x=reorder(ag_species2, nn), y=nn), fill=agcolor) + 
    # theme(axis.text.x = element_text(angle= 90)) + 
    geom_text(mapping = aes(x=ag_species2, y=nn+15, label = nn)) + 
    coord_flip() + 
    labs(x='  ', y='# of structures') +
    geom_text(mapping = aes(x=3, y = 200, label = total_label))
  
  outpdf(infile, 'agspecies_distribution_complex', width = 10, height = 15)
  
}

# run stuff
species_distribution_antibody()
# species_distribution_antigen()
# species_distribution_antigen_complex()
# species_distribution_antibody_complex()
# add complexes instead of paratopes 
