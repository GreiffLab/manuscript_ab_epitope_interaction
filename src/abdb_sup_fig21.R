# network aa physchem
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
library('wesanderson')
library('Peptides')
library('ggfortify')

names(wes_palettes)
theme_set(themeakbar())
my_spectral <- colorRampPalette(brewer.pal(8,'Spectral'))(14)
my_darjee = colorRampPalette(wes_palette('Darjeeling1'))(8)
print(my_darjee)
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
ragcolor_rgb = c(col2rgb('lightsalmon'))/255
ragcolor = rgb(ragcolor_rgb[1], ragcolor_rgb[2], ragcolor_rgb[3])
rabcolor_rgb = c(col2rgb('lightpink'))/255
rabcolor = rgb(rabcolor_rgb[1], rabcolor_rgb[2], rabcolor_rgb[3])

ppi_mcolor_rgb = c(col2rgb('orange'))/255
ppi_mcolor = rgb(ppi_mcolor_rgb[1], ppi_mcolor_rgb[2], ppi_mcolor_rgb[3])
ppi_pcolor_rgb = c(col2rgb('darkorange'))/255
ppi_pcolor = rgb(ppi_pcolor_rgb[1], ppi_pcolor_rgb[2], ppi_pcolor_rgb[3])

my_ccols = c('salmon',
             'deepskyblue1', 'deepskyblue2', 'deepskyblue3', 'deepskyblue4',
             'orange1', 'orange2', 'orange3', 'orange4')


my_ccols2 = c(abcolor3, 'lightpink',
              abcolor3, 'lightpink',
              abcolor3, 'lightpink',
              abcolor3, 'lightpink')

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


deltag_bar = function(){
  infile = 'abdb_outfiles_2019/foldx_deltag_split_chain.csv'
  df = drop_na(read_csv(infile))
  print(df)
  # df = df[df$deltag > 10,]
  # df$deltag = abs(df$deltag)
  print(df)
  mdf = df%>%group_by(abchain, segment,abresnum)%>% summarise(mean = median(deltag), sd = sd(deltag) )
  print(mdf)
  # stop()
  # xlabs = c(33,55,98,30,53,93,15,45,75,106,12,40,70,105)
  # xlabs = c(33,42,98,30,53,93,15,40,75,106,12,70,105)
  # stop()
  ggplot(data=mdf, mapping = aes(x=abresnum, y=mean, ymax=mean+(2*sd))) +
    geom_bar(stat = 'identity', mapping = aes(fill=segment)) +
    # geom_errorbar(mapping = aes(ymin=mean), color='grey') +
    facet_wrap(~ abchain, nrow = 2, scales = 'free') +
    scale_fill_manual(values = my_spectral) 
  outpdf(infile, 'bar')
}


deltag_hist = function(){
  infile = 'abdb_outfiles_2019/foldx_deltag_split_chain.csv'
  df = drop_na(read_csv(infile))
  print(df)
  df = df[df$deltag > 7,]
  print(df)
  countdf = count(df, segment, abchain)%>%group_by(abchain)%>% mutate(percent = round(n/sum(n)*100, 1))
  print(countdf)
  # stop()
  xlabs = c(33,55,98,30,53,93,15,45,75,106,12,40,70,105)
  ggplot(data=df, mapping = aes(x=deltag)) + 
    geom_histogram() + 
    facet_wrap(~ segment, scales = 'free')
  
  outpdf(infile, 'hist')
}


beat_deltag_bar = function(){
  infile = 'abdb_outfiles_2019/beatmusic_deltag_split_chain.csv'
  df = drop_na(read_csv(infile))
  print(df)
  ggplot(data=df, mapping = aes(x=abresnum, y=deltag)) +
    geom_bar(stat = 'identity', mapping = aes(fill=segment)) +
    # geom_errorbar(mapping = aes(ymin=mean), color='grey') +
    facet_wrap(~ abchain, nrow = 2, scales = 'free') +
    scale_fill_manual(values = my_spectral) 
  outpdf(infile, 'bar')
}

repair_unrepair = function(){
  infile = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin_foldx/try_repair/1A14_1_AS.fxout'
  df = read_delim(infile, delim = ' ')
  infile2 = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin_foldx/try_repair/1A14_1_Repair_AS.fxout'
  df2 = read_delim(infile2, delim = ' ')
  print(df)
  colnames(df) = c('res', 'pos', 'str1', 'str2', 'str3', 'str4','str5','deltag')
  print(df)
  print(df2)
  colnames(df2) = c('res', 'pos', 'str1', 'str2', 'str3', 'str4','str5','deltag')
  df$deltag2 = df2$deltag
  print(df)
  ggplot(data=df) + 
    geom_point(mapping = aes(x=deltag, y=deltag2)) 
  outpdf(infile, 'scatter')
}

# run stuff
# deltag_bar()
# deltag_hist()
# beat_deltag_bar()
# repair_unrepair()