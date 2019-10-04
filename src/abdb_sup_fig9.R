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



# run stuff
full_length_density = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx__full_length.csv'
  df = read_csv(infile)
  df[df=='H'] = 'Heavy'
  df[df=='L'] = 'Light'
  print(df)
  meds = df %>% group_by(abchain) %>% summarise(med_per_chain = median(motif_per_chain))
  maxs = df %>% group_by(abchain) %>% summarise(max_per_chain = max(motif_per_chain))
  med_labels = sprintf('Median: %s', meds$med_per_chain)
  print(meds)
  ggplot(data= df) + 
    # geom_density(mapping = aes(x=motif_per_chain, fill=abchain), color= NA) + 
    geom_bar(mapping = aes(x=motif_per_chain, fill=abchain), color= NA, position = 'dodge') + 
    scale_fill_manual(values = c(alpha(abcolor, 0.8), alpha(abcolor3,0.8))) + 
    # scale_color_manual(values = c(alpha(abcolor, 0.9), alpha(abcolor3,0.9))) + 
    # scale_color_manual(values = c(alpha(abcolor2, 0.6), alpha(abcolor3,0.9))) + 
    geom_text(data = meds, mapping = aes(x=c(1,1), y=c(300, 290), label= med_labels), size=5, color = c(abcolor, abcolor3)) + 
    labs(x= '# of interaction motifs per antibody', y= 'Density', fill = 'Antibody chain', color = "Antibody chain") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 25),
          legend.title = element_blank())
    # facet_wrap(~abchain)
  outpdf(infile, 'motif_perchain_density', height = 12, width = 7)
}


plot_pattern_bar_pie_full_length <- function(infile){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx__full_length.csv'
  df  =  read_csv(infile)
  print(df)
  countdf = count(df, abchain, gapstrstatus)
  countdf[countdf == 'H'] = 'Heavy'
  countdf[countdf == 'L'] = 'Light'
  g <- ggplot(data=countdf)
  g <- g + geom_bar(stat = 'identity', mapping = aes(x='', y=n, fill = gapstrstatus, color = gapstrstatus))
  g <- g + facet_wrap(~ abchain)
  g <- g + coord_polar(theta='y')
  g <- g + theme(axis.text.y = element_text(size = 1), axis.text.x = element_text(size = 8)) +
    labs(x=' ', y='# of continuous and discontinous motifs') +
    scale_fill_manual(values=c('red', 'black')) +
    scale_color_manual(values=c('red', 'black')) +
    labs (fill= 'Gap Status', color= 'Gap Status') + 
    theme(legend.position = 'none',
          axis.title = element_text(size= 23),
          strip.text = element_text(size = 20), 
          axis.text.x = element_text(size = 15)) +
    scale_y_continuous( breaks = c(150, 300, 450,600,750)) + 
    expand_limits(y=20 ) +
    # geom_text(mapping = aes(x=1, y= cumsum(raw_count), fill = gapstrstatus, color = gapstrstatus, label = df$label))
  outpdf(infile, 'pie', width = 8, height = 5)

}


plot_pattern_bar_radial_full_length <- function(infile){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx__full_length.csv'
  df  =  read_csv(infile)
  df[df == 'H'] = 'Heavy'
  df[df == 'L'] = 'Light'
  print(df)
  countdf = count(df, abchain, gapstrstatus)
  g <- ggplot(data=df)
  g <- g + geom_bar(mapping = aes(x=ab_motif, fill = gapstrstatus, color = gapstrstatus))
  g <- g + facet_wrap(~ abchain)
  g <- g + coord_polar(theta='y')
  g <- g + theme(axis.text.y = element_text(size = 1), axis.text.x = element_text(size = 8)) +
    labs(x=' ', y='# of shared continuous and discontinous motifs') +
    scale_fill_manual(values=c('red', 'black')) +
    scale_color_manual(values=c('red', 'black')) +
    labs (fill= 'Gap Status', color= 'Gap Status') + 
    theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + 
    theme(axis.title = element_text(size= 23),
          strip.text = element_text(size = 20), 
          axis.text.x = element_text(size = 15),
          legend.text = element_text(size=20),
          legend.title = element_blank()) +
    expand_limits(y=20 ) +
  outpdf(infile, 'radial', width = 8, height = 5)

}

# run stuff
full_length_density()
# plot_pattern_bar_pie_full_length()
# plot_pattern_bar_radial_full_length()


