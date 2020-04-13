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


plot.pattern.bar.pie <- function(infile){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv'
  df  =  read_csv(infile)
  df$label = scales::percent(df$normalized_count)
  print(df)
  g <- ggplot(data=df)
  g <- g + geom_bar(stat = 'identity', mapping = aes(x='', y= raw_count, fill = gapstrstatus, color = gapstrstatus))
  g <- g + facet_wrap(~ segment, ncol = 3)
  g <- g + coord_polar(theta='y')
  g <- g + theme(axis.text.y = element_text(size = 1), axis.text.x = element_text(size = 8)) +
    labs(x=' ', y='Number of continuous and discontinous motifs') +
    scale_fill_manual(values=c('red', 'black')) +
    scale_color_manual(values=c('red', 'black')) +
    labs (fill= 'Gap Status', color= 'Gap Status') + 
    theme(legend.position = 'none',
          axis.title = element_text(size = 23),
          strip.text = element_text(size = 20)) + 
    expand_limits(y=20 )
    # geom_text(mapping = aes(x=1, y= cumsum(raw_count), fill = gapstrstatus, color = gapstrstatus, label = df$label))
  outpdf(infile, 'pie_sup7', width = 8, height = 15)
}


plot.pattern.bar.radial <- function(infile){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv'
  # outname = strsplit(infile, split ='/')[[1]][2]
  # outname = paste0('abdb_figures_2019/',substr(outname, 1,nchar(outname)-4), '_radial.pdf') 
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
        axis.title = element_text(size = 23),
        strip.text = element_text(size = 20)) 
  outpdf(infile, 'radial_sup7', width = 8, height = 15)
}


# run stuff
# plot.pattern.bar.pie('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv')
# plot.pattern.bar.pie('abdb_outfiles_2019/respairs_absort_cutoff4_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv')
# plot.pattern.bar.pie('abdb_outfiles_2019/respairs_absort_cutoff6_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv')
plot.pattern.bar.radial('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv')
plot.pattern.bar.radial('abdb_outfiles_2019/respairs_absort_cutoff6_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv')
plot.pattern.bar.radial('abdb_outfiles_2019/respairs_absort_cutoff4_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv')