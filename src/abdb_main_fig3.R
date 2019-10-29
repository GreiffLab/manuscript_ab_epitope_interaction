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

pattern.bar.pie <- function(infile){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv'
  content <- read.csv(infile)
  print(summary(content))
  g <- ggplot(data=content)
  g <- g + geom_bar(stat = 'identity', mapping = aes(x='', y= raw_count, fill = gapstrstatus, color = gapstrstatus))
  g <- g + facet_wrap(~ segment, ncol = 3)
  g <- g + coord_polar('y', start=0)
  g <- g + theme(axis.text.y = element_blank(), axis.text.x = element_text(size=8)) +
    labs(x=' ', y='# of continuous and discontinous motifs') +
    scale_fill_manual(values=c('red', 'black')) +
    scale_color_manual(values=c('red', 'black')) +
    labs (fill= 'Gap Status', color= 'Gap Status') + 
    theme(legend.position = 'none',
          axis.title = element_text(size = 27),
          strip.text = element_text(size = 20))
  outname <- 'abdb_figures/shared_motifs_pie.pdf'
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_pie.pdf', outfigdir,inname)
  ggsave(outname, width=8, height=15)
  system(sprintf('open %s', outname))
}


pattern.bar.radial <- function(infile){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv'
  outname = strsplit(infile, split ='/')[[1]][2]
  outname = paste0('abdb_figures_2019/',substr(outname, 1,nchar(outname)-4), '_radial.pdf') 
  content <- read.csv(infile)
  print(summary(content))
  ggplot(data=content) + 
    geom_bar(stat = 'identity', mapping = aes(x=gap_patterns, y= raw_count, fill = gapstrstatus, color =
                                                gapstrstatus)) + 
    facet_wrap(~ segment, ncol = 3) + 
    coord_polar('y', start=0) +
    theme(axis.text = element_text(size = 10), axis.text.y = element_blank(), axis.ticks = element_blank()) +
    scale_fill_manual(values=c('red', 'black')) +
    scale_color_manual(values=c('red', 'black')) +
    # theme(panel.grid.major = element_line(colour = 'grey')) +
    labs (x= 'Motif', y='Number of shared motifs', fill= 'Gap Status', color= 'Gap Status') +
    theme(legend.position = 'none',
          axis.title = element_text(size = 27),
          strip.text = element_text(size = 20))
  ggsave(outname, width=8, height=15)
  system(sprintf('open %s', outname))
}

shared_structures = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil_pc_reglen.csv'
  df = read_csv(infile)
  countdf  = count(df, pdbid, ab_motif)
  countdf = count(countdf, ab_motif)
  countdf2  = count(df, pdbid, ag_motif)
  countdf2 = count(countdf2, ag_motif)
  print(countdf)
  outdf = data_frame()
  for (i in seq(2,20)){
    print(i)
    sdf = countdf[countdf$n>=i,]
    # print(dim(sdf))
    newrow = data_frame(nmotif = dim(sdf)[1], cutoff = i, source = 'Paratope')
    outdf = rbind(outdf, newrow)
    sdf2 = countdf2[countdf2$n>i,]
    # print(dim(sdf))
    newrow2 = data_frame(nmotif = dim(sdf2)[1], cutoff = i, source = 'Epitope')
    outdf = rbind(outdf, newrow2)
    # stop()
  }
  outdf$source = factor(outdf$source, level = c('Paratope', 'Epitope'))
  ggplot(data=outdf) +
    geom_bar(mapping = aes(x=cutoff, y=nmotif, fill=source), stat = 'identity') + 
    facet_wrap(~source, scales = 'free')+ 
    geom_text(data=outdf, mapping = aes(x=cutoff, y=nmotif*1.1, label = nmotif)) + 
    labs(x= 'Cut-offs for number of structures \n(2-20 structures)', y = 'Motifs found in a cut-off') + 
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 25),
          legend.title = element_blank(), 
          # legend.text = element_text(size=25),
          legend.position = 0,
          strip.text = element_text(size = 25),) +
    scale_fill_manual(values = c(abcolor3, agcolor))
  outpdf(infile, 'shared_structures', width = 20)
}

# run stuff
# pattern.bar.pie('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv')
# pattern.bar.pie('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_count.csv')
# pattern.bar.radial('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv')
# pattern.bar.radial('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_count.csv')
shared_structures()
