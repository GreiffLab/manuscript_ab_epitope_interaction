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
ragcolor_rgb = c(col2rgb('lightsalmon'))/255
ragcolor = rgb(ragcolor_rgb[1], ragcolor_rgb[2], ragcolor_rgb[3])
rabcolor_rgb = c(col2rgb('lightpink'))/255
rabcolor = rgb(rabcolor_rgb[1], rabcolor_rgb[2], rabcolor_rgb[3])


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


ppi_mcolor_rgb = c(col2rgb('orange'))/255
ppi_mcolor = rgb(ppi_mcolor_rgb[1], ppi_mcolor_rgb[2], ppi_mcolor_rgb[3])
ppi_pcolor_rgb = c(col2rgb('darkorange'))/255
ppi_pcolor = rgb(ppi_pcolor_rgb[1], ppi_pcolor_rgb[2], ppi_pcolor_rgb[3])

my_ccols = c('salmon',
             'deepskyblue1', 'deepskyblue2', 'deepskyblue3', 'deepskyblue4',
             'orange1', 'orange2', 'orange3', 'orange4')


my_ccols2 = c('salmon',
             abcolor3, abcolor3, abcolor3, abcolor3,
             ppi_mcolor, ppi_mcolor, ppi_mcolor, ppi_mcolor)

my_ccols3 = c(ggplot2::alpha('lightpink', 1), ggplot2::alpha(abcolor,0.5),
              ggplot2::alpha('lightpink', 1), ggplot2::alpha(abcolor,0.5), 
              ggplot2::alpha('lightpink',1), ggplot2::alpha(abcolor,0.5), 
              ggplot2::alpha('red', 0.5), ggplot2::alpha(abcolor3, 1))
print(my_ccols3)
# stop()

dl_ld = function(){
  infile = 'abdb_outfiles_2019/merged_eval_files.csv'
  # infile = 'abdb_outfiles_2019/merged_eval_files10k.csv'
  df = read_csv(infile)
  dfmed = df %>% group_by(rep, use_case, data_tag, exp_tag) %>% summarise(med = median(ldnorm), sd = sd(ldnorm))
  dfmed = dfmed %>% group_by(use_case, data_tag,exp_tag) %>% summarise(me = mean(med), se = sd(med)/sqrt(length(med)), len = length(med))
  print(dfmed[order(dfmed$me),])
  print(dfmed)
  # stop()
  xticklabels=  c('Motif E --> P', 'Motif + position E --> P',
                 'Motif P --> E', 'Motif + position P --> E',
                 'Sequence E --> P', 'Sequence P --> E')
  xtickbreaks = unique(dfmed$use_case)
  print(xticklabels)
  print(xtickbreaks)
  # stop()
  ggplot(data=dfmed, mapping = aes(x=use_case, y=me, ymin=me, ymax=me+2*se, fill = exp_tag)) + 
    geom_errorbar(mapping = aes(ymin=me+2*se), position = position_dodge(0.9), width = 0.1) +
    geom_linerange(position = position_dodge(0.9)) +
    geom_bar( position = position_dodge(), stat = 'identity') + 
    geom_text(aes(y= me-0.02, label = round(me, 2)), position = position_dodge(0.9)) + 
    labs(y= 'Median prediction error', x= 'Use case', fill = 'Experiment type') + 
    # scale_fill_discrete(labels =c('Randomized pairs', 'Observed pairs')) + 
    scale_fill_manual(values = c(alpha('lightsalmon', 0.5), 'lightsalmon'), labels =c('Randomized pairs', 'Observed pairs')) +
    scale_x_discrete(breaks = xtickbreaks, labels = xticklabels)
  outpdf(infile, 'bar', height = 7, width = 12)
}



plot_sl_dl_summary = function(){
  # note to self: the merge sl dl abdb_prepdata_sup_fig18, now also in prepdata_sup_fig5
  infile = 'abdb_outfiles_2019/sl_dl_evalsummary.csv'
  df = read_csv(infile)
  df = df[df$source == 'abdb',]
  # print(unique(df$exp_tag))
  print(df$repldexactmea)
  # testv = 2*df$repldexactmea
  # print(testv)
  df$exp_tag = factor(df$exp_tag, levels = c('control_marginal_proba', 'exp_marginal_proba',
                                              'control_cond_proba', 'exp_cond_proba',
                                              'control_cond_proba_with_prior', 'exp_cond_proba_with_prior',
                                              'control_NMT', 'exp_NMT'))
  print('yo')
  print(df$use_case)
  df$use_case = factor(df$use_case, levels = c('motif_paraepi', 'motif_paraepipos', 'motif_epipara', 'motif_epiparapos',
                                               'seq_paraepi', 'agg_paraepi', 'seq_epipara', 'agg_epipara'))
  print(df$exp_tag)
  intags = c('control_cond_proba_with_prior', 'exp_cond_proba_with_prior', 'exp_NMT', 'control_NMT')
  df = df %>% filter(exp_tag %in% intags)
  print(dim(df))
  # stop()
  xticklabels=  c('Motif P to E', 'Motif + pos P to E',
                 'Motif E to P', 'Motif + pos E to P',
                 'Seq P to E', 'Agg P to E',
                 'Seq E to P', 'Agg E to P')
  # xtickbreaks = unique(dfmed$use_case)
  print(xticklabels)
  barlabels = c()
  for (tag in df$exp_tag){
    if (grepl('NMT', tag)){
      tag = 'DEEP' 
      barlabels = c(barlabels, tag)
    }
    if (grepl('cond_proba_with', tag)){
      tag = 'SHALLOW' 
      barlabels = c(barlabels, tag)
    }
    if (grepl('cond_proba', tag)){
      tag = 'Shallow 2' 
      barlabels = c(barlabels, tag)
    }
    if (grepl('marginal', tag)){
      tag = 'Shallow 1' 
      barlabels = c(barlabels, tag)
    }
  }
  print(barlabels)
  df$exp_tag2 = barlabels
  print(df)
  # stop()
  # print(xtickbreaks)
  # print(df)
  ggplot(data=df, mapping = aes(x=use_case, y=repldnormmea, ymin=repldnormmea, ymax=repldnormmea+2*repldnormse, fill = exp_tag)) + 
    geom_errorbar(mapping = aes(ymin=repldnormmea+2*repldnormse), color = 'gray', position = position_dodge(0.9), width = 0.1) +
    geom_linerange(position = position_dodge(0.9), color = 'gray') +
    geom_bar( position = position_dodge(), stat = 'identity') + 
    geom_text(aes(y= repldnormmea-0.01, label = round(repldnormmea, 2)), size= 4, color = 'black', position = position_dodge(0.9)) + 
    geom_text(aes(y=0.005, label = exp_tag2), angle=90, size= 5, color='black', position = position_dodge(0.9), hjust=0) + 
    labs(y= 'Mean prediction error', x= 'Use case', fill = 'Experiment type') + 
    scale_fill_manual(values=my_ccols3) +
    theme(legend.position = 0) + 
    scale_x_discrete(labels = xticklabels) +
    theme(axis.title = element_text(size = 30, color = 'gray42'), 
          axis.text = element_text(size=22, color='gray42'))
  outpdf(infile, 'norm_bar', height = 11, width = 25)
  stop()
  ggplot(data=df, mapping = aes(x=use_case, y=repldexactmea, ymin=repldexactmea, ymax=repldexactmea+testv, fill = exp_tag)) + 
    geom_errorbar(mapping = aes(ymin=repldexactmea+testv), color= 'gray',position = position_dodge(0.9), width = 0.1) +
    geom_linerange(position = position_dodge(0.9), color = 'gray') +
    geom_bar( position = position_dodge(), stat = 'identity') + 
    geom_text(aes(y= repldexactmea-0.01, label = round(repldexactmea, 2)), color = 'white', size= 3, position = position_dodge(0.9)) + 
    geom_text(aes(y=0.005, label = exp_tag), size= 10, color='white', angle = 90, hjust=0, position = position_dodge(0.9)) + 
    labs(y= 'Mean prediction error', x= 'Use case', fill = 'Experiment type') + 
    scale_fill_manual(values=my_ccols3) + 
    theme(legend.position = 0) + 
    # scale_x_discrete(labels = xticklabels) + 
    theme(axis.title = element_text(size = 20, color = 'gray42'), axis.text = element_text(size=15, color='gray42'))
  outpdf(infile, 'exact_bar', height = 7, width = 20)
  
}


ld_bin = function(x){
  x[x>0] = 1
  return = x
}

dl_ld_binary = function(){
  # output ldexact and ldnorm centers median and mean.
  infile = 'abdb_outfiles_2019/merged_eval_files.csv'
  # infile = 'abdb_outfiles_2019/merged_eval_files10k.csv'
  df = read_csv(infile)
  print(df)
  bins = ld_bin(df$ldnorm) # binary ld
  print(sum(bins))
  print(dim(df))
  df$ldexact = bins
  dfmed = df %>% group_by(rep, use_case, data_tag, exp_tag) %>% summarise(ldnorm_mea = mean(ldnorm), 
                                                                          ldnorm_sd = sd(ldnorm),
                                                                          ldexact_mea = mean(ldexact),
                                                                          ldexact_sd = sd(ldexact))
  print(dfmed)
  print(dfmed[order(dfmed$ldnorm_mea),])
  dfmed = dfmed %>% group_by(use_case, data_tag,exp_tag) %>% summarise(repldnormmea = mean(ldnorm_mea), 
                                                                       repldnormse = sd(ldnorm_mea)/sqrt(length(ldnorm_mea)), 
                                                                       ldnormreps = length(ldnorm_mea),
                                                                       repldexactmea = mean(ldexact_mea), 
                                                                       repldexactse = sd(ldexact_mea)/sqrt(length(ldexact_mea)), 
                                                                       ldexactreps = length(ldexact_mea))
  xticklabels=  c('Motif E --> P', 'Motif + position E --> P',
                 'Motif P --> E', 'Motif + position P --> E',
                 'Sequence E --> P', 'Sequence P --> E')
  xtickbreaks = unique(dfmed$use_case)
  print(xticklabels)
  print(xtickbreaks)
  dfmed[dfmed=='control'] = 'control_NMT'
  dfmed[dfmed=='exp'] = 'exp_NMT'
  dfmed['source'] = rep('abdb', dim(dfmed)[1])
  print(dfmed)
  # stop()
  write_csv(dfmed, 'abdb_outfiles_2019/eval_summary.csv')
  stop()
}
  
agg_dl = function(){
  infile = 'abdb_outfiles_2019/agg_merged_eval_files.csv'
  df = read_csv(infile)
  print(df)
  bins = ld_bin(df$ldnorm) # binary ld
  print(sum(bins))
  print(dim(df))
  df$ldexact = bins
  dfmed = df %>% group_by(rep, use_case, data_tag, exp_tag) %>% summarise(ldnorm_mea = mean(ldnorm), 
                                                                          ldnorm_sd = sd(ldnorm),
                                                                          ldexact_mea = mean(ldexact),
                                                                          ldexact_sd = sd(ldexact))
  print(dfmed)
  print(dfmed[order(dfmed$ldnorm_mea),])
  dfmed = dfmed %>% group_by(use_case, data_tag,exp_tag) %>% summarise(repldnormmea = mean(ldnorm_mea), 
                                                                       repldnormse = sd(ldnorm_mea)/sqrt(length(ldnorm_mea)), 
                                                                       ldnormreps = length(ldnorm_mea),
                                                                       repldexactmea = mean(ldexact_mea), 
                                                                       repldexactse = sd(ldexact_mea)/sqrt(length(ldexact_mea)), 
                                                                       ldexactreps = length(ldexact_mea))
  
  # print(dfmed[order(dfmed$repldnormmea),])
  dfmed[dfmed=='control'] = 'control_NMT'
  dfmed[dfmed=='exp'] = 'exp_NMT'
  dfmed[dfmed=='seq'] = 'agg'
  dfmed[dfmed=='seq_epipara'] = 'agg_epipara'
  dfmed[dfmed=='seq_paraepi'] = 'agg_paraepi'
  dfmed['source'] = rep('abdb', dim(dfmed)[1])
  print(dfmed)
  # stop()
  write_csv(dfmed, 'abdb_outfiles_2019/agg_eval_summary.csv')
}

# run stuff
# dl_ld()
# dl_ld_binary()
# agg_dl()
plot_sl_dl_summary()
# 