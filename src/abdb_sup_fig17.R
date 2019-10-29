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



dl_length_correlation = function(infile) {
  # infile = 'abdb_outfiles_2019/merged_eval_files.csv'
  # infile = 'abdb_outfiles_2019/merged_eval_files_paraepi.csv'
  # infile = 'abdb_outfiles_2019/merged_eval_files10k.csv'
  df = read_csv(infile)
  ### output control and exp
  # dfpe = df[(df$use_case == 'motif_paraepi') & (df$rep == 10) & (df$exp_tag == 'control'),]
  # print(dfpe)
  # write_csv(dfpe, 'abdb_outfiles_2019/merged_eval_files_paraepi_control.csv')
  # dfpe = df[(df$use_case == 'motif_paraepi') & (df$rep == 10) & (df$exp_tag == 'exp'),]
  # print(dfpe)
  # write_csv(dfpe, 'abdb_outfiles_2019/merged_eval_files_paraepi_exp.csv')
  # stop()
  ## end output control and exp
  print(df)
  ggplot(data = df, mapping = aes(x=lenpredpara, y=lenpara)) + 
    geom_point(color = alpha('lightsalmon', 0.1), size = 2) +
    geom_density2d(color = alpha('lightsalmon', 0.9)) + 
    facet_grid(embeding_dim ~ units) + 
    geom_smooth(method = 'lm', color = alpha('lightslateblue', 0.5)) + 
    labs(x= 'Length predicted output', y = 'Length observed output') + 
    theme(axis.title = element_text(size = 40),
          strip.text = element_text(size = 28))
  outpng(infile, 'length_correlation', height = 15, width = 18)
}


dl_length_correlation_bar = function(){
  # infile = 'abdb_outfiles_2019/merged_eval_files_paraepi.csv'
  # infile = 'abdb_outfiles_2019/merged_eval_files10k.csv'
  infile = 'abdb_outfiles_2019/merged_eval_files.csv'
  infile2 = 'abdb_outfiles_2019/agg_merged_eval_files.csv'
  # infile = 'abdb_outfiles_2019/merged_eval_files10k.csv'
  df = read_csv(infile)
  df2 = read_csv(infile2)
  df2[df2 == 'seq_epipara'] = 'agg_epipara'
  df2[df2 == 'seq_paraepi'] = 'agg_paraepi'
  print(unique(df2$use_case))
  df = rbind(df, df2)
  print(df)
  # stop()
  df$use_case = factor(df$use_case, levels = c('motif_paraepi', 'motif_paraepipos', 'motif_epipara', 'motif_epiparapos',
                                               'seq_paraepi', 'agg_paraepi', 'seq_epipara', 'agg_epipara'))
  print(df)
  # stop()
  xticklabels=  c('Motif P to E', 'Motif + pos P to E',
                 'Motif E to P', 'Motif + pos E to P',
                 'Seq P to E', 'Agg P to E',
                 'Seq E to P', 'Agg E to P')
  # stop()
  dfmed = df %>% group_by(rep, use_case, data_tag, exp_tag, embeding_dim, units) %>% summarise(med = cor(lenpredpara, lenpara), sd = sd(ldnorm))
  dfmed = dfmed %>% group_by(use_case, data_tag,exp_tag, embeding_dim, units) %>% summarise(me = mean(med), se = sd(med)/sqrt(length(med)), len = length(med))
  print(dfmed[order(-dfmed$me),])
  dfmed = dfmed %>% group_by(use_case, data_tag,exp_tag) %>% summarise(ma = max(me, na.rm = TRUE), se = sd(me, na.rm = TRUE)/sqrt(length(me)), len = length(me))
  print(dfmed[order(-dfmed$ma),])
  print(dfmed)
  is.na(dfmed) <- sapply(dfmed, is.infinite)
  dfmed[is.na(dfmed)] = 0 
  print(dfmed)
  # dfmed[dfmed==NA] = 0
  # stop()
  # xticklabels=  c('Motif E --> P', 'Motif + position E --> P',
  #                'Motif P --> E', 'Motif + position P --> E',
  #                'Sequence E --> P', 'Sequence P --> E')
  # xtickbreaks = unique(dfmed$use_case)
  # stop()
  ggplot(data=dfmed, mapping = aes(x=use_case, y=ma, ymin=ma, ymax=ma+2*se, fill = exp_tag)) + 
    geom_bar( position = position_dodge(), stat = 'identity') +
    geom_errorbar(mapping = aes(ymin=ma+2*se), position = position_dodge(0.9), width = 0.1) +
    geom_linerange(position = position_dodge(0.9)) +
    geom_text(aes(y= ma-0.02, label = round(ma, 2)), position = position_dodge(0.9)) +
    labs(y= 'Predicted vs observed length \n(Pearson cor.)', x= 'Use case', fill = 'Experiment type') + 
    # scale_fill_discrete(labels =c('Randomized pairs', 'Observed pairs')) + 
    scale_fill_manual(values = c(alpha('lightsalmon', 0.5), 'lightsalmon'), labels =c('Randomized pairs', 'Observed pairs')) +
    # scale_x_discrete(breaks = xtickbreaks, labels = xticklabels) + 
    scale_x_discrete(labels = xticklabels) + 
    theme(axis.title = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          legend.text = element_text(size = 20),
          legend.title = element_blank())
  outpdf(infile, 'corr_bar', height = 7, width = 15)
}


## run stuff
# dl_length_correlation('abdb_outfiles_2019/merged_eval_files_paraepi_exp.csv')
# dl_length_correlation('abdb_outfiles_2019/merged_eval_files_paraepi_control.csv')
dl_length_correlation_bar()



