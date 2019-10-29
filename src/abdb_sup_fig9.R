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
  print(df)
  # stop()
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
    geom_text(data = meds, mapping = aes(x=c(1,1), y=c(300, 290), label= med_labels), size=7, color = c(abcolor, abcolor3), hjust=0) + 
    labs(x= '# of interaction motifs per antibody', y= '# of antibody-antigen complexes', fill = 'Antibody chain', color = "Antibody chain") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 25),
          legend.title = element_blank())
    # facet_wrap(~abchain)
  outpdf(infile, 'motif_perchain_density', height = 12, width = 7)
}



full_length_residues = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx__full_length.csv'
  df = read_csv(infile)
  # stop()
  df[df=='H'] = 'Heavy'
  df[df=='L'] = 'Light'
  print(df)
  df = df %>% group_by(pdbid, abchain) %>% mutate(res_len = nchar(paratope))
  meds = df %>% group_by(abchain) %>% summarise(med_per_chain = median(res_len))
  maxs = df %>% group_by(abchain) %>% summarise(max_per_chain = max(motif_per_chain))
  med_labels = sprintf('Median: %s', meds$med_per_chain)
  print(meds)
  # stop()
  ggplot(data= df) + 
    # geom_density(mapping = aes(x=motif_per_chain, fill=abchain), color= NA) + 
    geom_bar(mapping = aes(res_len, fill=abchain), color= NA, position = 'dodge') + 
    scale_fill_manual(values = c(alpha(abcolor, 0.8), alpha(abcolor3,0.8))) + 
    # scale_color_manual(values = c(alpha(abcolor, 0.9), alpha(abcolor3,0.9))) + 
    # scale_color_manual(values = c(alpha(abcolor2, 0.6), alpha(abcolor3,0.9))) + 
    geom_text(data = meds, mapping = aes(x=c(1,1), y=c(120, 116), label= med_labels), size=7, color = c(abcolor, abcolor3), hjust=0) + 
    labs(x= '# of interacting residues per antibody', y= '# of antibody-antigen complexes', fill = 'Antibody chain', color = "Antibody chain") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 25),
          legend.title = element_blank())
    # facet_wrap(~abchain)
  outpdf(infile, 'residues_perchain_density', height = 12, width = 7.5)
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


ab_ag_motiflen_correlation  = function(){
  # examine ab ag angle per region
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged.csv'
  df = read.csv(infile)
  print(head(df))
  df = drop_na(df)
  segments = sort(unique(df$segment))
  print(segments[1])
  cors_data = sapply(1:length(segments), function(x) cor.test(subset(df, df$segment == segments[x])$ab_motiflen,subset(df, df$segment == segments[x])$ag_motiflen))
  cors = sapply(1:dim(cors_data)[2], function(x) round(cors_data[,x]$estimate,2))
  print(cors_data[,1])
  pvals = sapply(1:dim(cors_data)[2], function(x) signif(cors_data[,x]$p.value, digits = 1))
  # cors_pearson = df %>% group_by(segment) %>% summarize(cor = cor(p_angle, e_angle, method='pearson', use = 'complete.obs'))
  cors_pearson = df %>% group_by(segment) %>% summarize(cor = cor(ab_motiflen, ag_motiflen, method='pearson', use = 'complete.obs'))
  cors_spearman = df %>% group_by(segment) %>% summarize(cor = cor(ab_motiflen, ag_motiflen, method='spearman', use = 'complete.obs'))
  print(cors_pearson)
  print(cors_spearman)
  label_pearson = paste0('Pearson Cor. ', round(cors_pearson$cor, 2))
  label_spearman = paste0('Spearman Cor. ', round(cors_spearman$cor, 2))
  print(label_pearson)
  ggplot(data=df, mapping = aes(x=ab_motiflen, y=ag_motiflen)) + 
    geom_point(mapping = aes(x = ab_motiflen, y= ag_motiflen), color = alpha(my_spectral[3],0.5), size=0.5) + 
    geom_text(data = cors_pearson, mapping = aes(x=2, y=20, label=label_pearson), hjust=0) +
    geom_text(data = cors_pearson, mapping = aes(x=2, y=17, label=label_spearman), hjust=0) +
    facet_wrap(~ segment, ncol = 3) + 
    geom_smooth(method = 'lm', mapping = aes(x=ab_motiflen, y= ag_motiflen), color=alpha(abcolor, 0.3)) + 
    geom_density2d(color = alpha(my_spectral[3],0.7)) +
    labs(x='Interaction motif length (paratope)', y='Interaction motif length (epitope)') + 
    theme(axis.title = element_text(size= 20),
          strip.text = element_text(size = 17), 
          axis.text.x = element_text(size = 15),
          legend.text = element_text(size=20),
          legend.title = element_blank()) 
  outname = 'abdb_figures/ab_ag_motiflen_correlation.pdf'
  ggsave(outname, height = 9, width=6)
  system(sprintf('open %s', outname))
  
}

region_len_corr = function(var1, var2){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil_pc_reglen.csv'
  df = read_csv(infile)
  print(df)
  if (var2 == 'plen'){ 
    ylabel = 'Paratope length'
    point_col = abcolor3}
  if (var2 == 'ab_motiflen'){ 
    ylabel = 'Paratope motif length' 
    point_col = abcolor3}
  if (var2 == 'epitope_len'){ 
    ylabel = 'Epitope length'
    point_col = agcolor}
  if (var2 == 'ag_motiflen'){ 
    ylabel = 'Epitope motif length'
    point_col = agcolor}
  xlabel = 'Region length'
  corr = cor(df$region_len, df$plen)
  print(corr)
  segments = unique(df$segment)
  print(segments)
  corrdf = data_frame()
  for (segment in segments){
    sdf = df[df$segment==segment,]
    cor  = cor(get(var1, sdf), get(var2, sdf))
    cor = round(cor, 2)
    if (is.na(cor)){cor = 'NA'}
    # cor = format(cor, scientific = TRUE, digits = 2)
    print(cor)
    cdf  = data_frame(segment= segment, cor = cor, x = min(get(var1,sdf)), y = max(get(var2, sdf)))
    corrdf = rbind(corrdf, cdf)
    # stop()
  }
  print(corrdf)
  # stop()
  ggplot(data = df) + 
    geom_point(mapping = aes(x=get(var1, df) , y=get(var2, df)), alpha = 0.2, color = point_col) + 
    facet_wrap(~segment, ncol=3, scales = 'free') + 
    scale_x_continuous(breaks = pretty_breaks(n=4)) + 
    labs(x= xlabel, y=ylabel) +
    geom_text(data = corrdf, mapping = aes(x=x-4, y=y*0.8, label = sprintf('Pear. cor.: %s', cor)), hjust=0,vjust=0, size = 4) +
    # ylim(c(-1,20)) + 
    theme(axis.title = element_text(size=20),
          strip.text = element_text(size =17))
  tag= sprintf('%s_%s_corr', var1, var2)
  print(tag)
  outpdf(infile, tag, width = 6)
}


# run stuff
# full_length_density()
# full_length_residues()
# plot_pattern_bar_pie_full_length()
# plot_pattern_bar_radial_full_length()
# ab_ag_motiflen_correlation()
region_len_corr('region_len', 'plen')
region_len_corr('region_len', 'epitope_len')
# region_len_corr('region_len', 'ab_motiflen')
# region_len_corr('region_len', 'ag_motiflen')

