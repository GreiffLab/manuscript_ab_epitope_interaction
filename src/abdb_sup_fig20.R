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


physchem_pca = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil_pc.csv'
  df = read_csv(infile)
  res = aaComp('BLADIHELL')
  print(res)
  stop()
  # df = df[(df$plen >1) & (df$epitope_len > 1),]
  # df$ab_motif[df$ab_motif != 'XXXX'] = '0'
  # df$ab_motif[df$ab_motif == 'XXXX'] = '1'
  print(df)
  df2 = rbind(df[df$ab_motif=='XXX',], df[df$ab_motif=='XXXX',], df[df$ab_motif=='XXXX1X',])
  # df2 = rbind(df[df$ab_motif=='XXX',], df[df$ab_motif=='XXXX',], df[df$ab_motif=='XXX1X',], df[df$ab_motif=='XXXXXX',],
  #             df[df$ab_motif=='XXXX1X',], df[df$ab_motif=='XXXXXXX',], df[df$ab_motif=='XX',], df[df$ab_motif=='XX1X',],
  #             df[df$ab_motif=='X1X',], df[df$ab_motif=='XXXXX',])
  df = df2
  # stop()
  dfpc = data_frame()
  for (i in 1:dim(df2)[1]){
    paratope = df$paratope[i]
    epitope = df$epitope[i]
    physchem = unlist(aaComp(epitope))[1:18]
    # physchem = unlist(aaComp(paratope))[10:18]
    # physchem2 = unlist(aaComp(epitope))[1:18]
    # physchem = c(physchem, physchem2)
    n = length(physchem)
    dfpc[i,1] = df$ab_motif[i]
    dfpc[i,2] = epitope
    dfpc[i, 3:(n+3-1)] = physchem
    # print(dfpc)
    # stop()
  }
  print(dfpc)
  # stop()
  pcs = prcomp(dfpc[,c(3:(n+3-1))], center = TRUE, scale. = TRUE)
  print(pcs)
  summary(pcs)
  autoplot(pcs, data=dfpc, fill = 'V1', 0.2, shape = 21, size=5, ) +
    theme(axis.title = element_text(size=10))
  outpdf(infile, 'pca', height = 7, width = 7)
}


# run stuff
physchem_pca()
