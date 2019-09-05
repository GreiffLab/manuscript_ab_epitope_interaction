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

abdb_bnaber = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber.csv'
  df = read_csv(infile)
  print(df)
  print(count(df, pdbid, bn_status) %>% count(bn_status))
  # stop()
  bndfh  = subset(df, (df$bn_status == 'BN') & (df$abchain == 'H'))
  bndfl  = subset(df, (df$bn_status == 'BN') & (df$abchain == 'L'))
  bnmotifh = bndfh$ab_motif
  bnmotifl = bndfl$ab_motif
  abmotif = df$ab_motif
  x = list(bnmotifh, bnmotifl, abmotif)
  partitions = get.venn.partitions(x)
  print(partitions$..values..)
  # stop()
  outname  = 'abdb_figures_2019/bnh_bnl_paratope_venn.tiff'
  venn.plot <- venn.diagram(x=x, filename=outname,
                            output = FALSE ,
                            imagetype="tiff" ,
                            # height = 800 , 
                            # width = 900 , 
                            resolution = 300,
                            compression = "lzw",
                            lwd = 2,
                            lty = 'blank',
                            fill = c('red', 'red', abcolor3),
                            alpha=c(0.5,0.5, 0.4),
                            cex = 5, 
                            cat.fontface=3,
                            cat.cex = 3,
                            cat.dist = c(-0.05,-0.05, -0.06),
                            category.names=c("bNAb heavy", "bNAb light", "Non-bNAb")
  ) 
  system(sprintf('open %s', outname))
  stop()
  ggplot(data= df ) + 
    geom_bar(mapping = aes(x=ab_motif)) + 
    facet_wrap(~ BN_status, ncol = 2, scales = 'free') + 
    coord_flip()
  
  outpdf(infile, 'bar', width = 14)
}

abdb_bnaber_top20 = function() {
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber.csv'
  df = read_csv(infile)
  df[df == 'BN'] = 'bNAb'
  df[df == 'Non-BN'] = 'Non-bNAb'
  n = 20
  countdf = count(df, ab_motif, bn_status, sort = TRUE)
  print(countdf)
  countdf2 = count(df, ab_motif, sort = TRUE)
  topn_abmotif = countdf2[1:n,]$ab_motif
  topndf = countdf[countdf$ab_motif %in% topn_abmotif,]
  print(unique(topndf$ab_motif))
  print(topn_abmotif)
  ggplot(data =topndf) + 
    geom_bar(mapping = aes(x=reorder(ab_motif, -n), y=n, fill=bn_status), stat='identity', position='dodge') + 
    geom_text(data = topndf, mapping = aes(x=ab_motif, y = n+20, label = topndf$n, fill=bn_status), position = position_dodge(width=1)) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    scale_fill_manual(values = c(abcolor3, alpha(abcolor3, 0.5))) + 
    labs(x = "Top 20 paratope motifs", y = '# of motifs', fill = 'Broadly neutralizing status')
  outpdf(infile, 'bar')
}

bnaber_abdb_diff = function() {
  infile = 'abdb_outfiles_2019/bnaber_abdb_diff.csv'
  df = read_csv(infile)
  countdf =  count(df, poly_count)
  print(countdf)
  print(df[df$poly_count > 2,])
  stop()
  ggplot(data = df) + 
    geom_bar(mapping = aes(x=poly_count))
  
}

# run stuff
abdb_bnaber()
# abdb_bnaber_top20()
# bnaber_abdb_diff()
