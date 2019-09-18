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

ppi_mcolor_rgb = c(col2rgb('orange'))/255
ppi_mcolor = rgb(ppi_mcolor_rgb[1], ppi_mcolor_rgb[2], ppi_mcolor_rgb[3])
ppi_pcolor_rgb = c(col2rgb('darkorange'))/255
ppi_pcolor = rgb(ppi_pcolor_rgb[1], ppi_pcolor_rgb[2], ppi_pcolor_rgb[3])


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

ld_bin = function(x){
  x[x>0] = 1
  return = x
}

dl_ld_binary = function(){
  infile = 'abdb_outfiles_2019/merged_eval_files.csv'
  # infile = 'abdb_outfiles_2019/merged_eval_files10k.csv'
  df = read_csv(infile)
  print(df)
  bins = ld_bin(df$ldnorm) # binary ld
  print(sum(bins))
  print(dim(df))
  df$ldnorm = bins
  dfmed = df %>% group_by(rep, use_case, data_tag, exp_tag) %>% summarise(med = mean(ldnorm), sd = sd(ldnorm))
  dfmed = dfmed %>% group_by(use_case, data_tag,exp_tag) %>% summarise(me = mean(med), se = sd(med)/sqrt(length(med)), len = length(med))
  print(dfmed[order(dfmed$me),])
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
    labs(y= 'Mean prediction error', x= 'Use case', fill = 'Experiment type') + 
    # scale_fill_discrete(labels =c('Randomized pairs', 'Observed pairs')) + 
    scale_fill_manual(values = c(alpha('lightsalmon', 0.5), 'lightsalmon'), labels =c('Randomized pairs', 'Observed pairs')) +
    scale_x_discrete(breaks = xtickbreaks, labels = xticklabels)
  outpdf(infile, 'absolute_bar', height = 7, width = 12)
}



chao_diversity = function(){
  #compute diversity estimate using chao2
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv'
  df = read_csv(infile)
  print(df)
  # stop()
  ## retrofit column to match ab-ag data
  df = df %>% rename(ab_motif = gap_pattern1,
                  ag_motif = gap_pattern2)
  ##
  dfab = df %>% group_by(ab_motif) %>%tally()
  dfab = dfab[order(-dfab$n),]
  # chao_ab = chao1(dfab$n)
  # print(chao_ab)
  # print(dfab)
  dfag = df %>% group_by(ag_motif) %>%tally()
  dfag = dfag[order(-dfag$n),]
  i.out = iNEXT(dfab$n, datatype = 'abundance')
  chao_abinext = ChaoRichness(dfab$n)
  chao_aginext = ChaoRichness(dfag$n)
  # print(chao_ab)
  # print(chao_abinext)
  # chao_ag = chao1(dfag$n)
  # print(chao_ag)
  divdf = data.frame('divs' = c(chao_abinext$Estimator, chao_aginext$Estimator, chao_abinext$Observed, chao_aginext$Observed), 'type' = c('Estimated (Chao)', 'Estimated (Chao)', 'Observed', 'Observed'), 'source' = c('Paratope', 'Epitope', 'Paratope', 'Epitope'))
  print(divdf)
  divdf$source = factor(c('PPI motif', 'PPI motif partner'), levels = c('PPI motif', 'PPI motif partner'))
  ytext = divdf$divs+(0.08*divdf$divs)
  ytext = round(ytext)
  ggplot(data = divdf) +
    geom_bar(stat = 'identity', mapping = aes(reorder(type,divs), divs, fill = source), position = 'dodge') + 
    scale_fill_manual(values = c(ppi_mcolor, ppi_pcolor), name = "Motif source", labels = c('PPI motif', 'PPI motif partner')) + 
    facet_wrap(~source, scales='free') + 
    geom_text(mapping = aes(x=type, y=ytext, label=round(divs))) +
    labs(x='', y= 'Interaction motif diversity')
  outpdf(infile, 'chao')
  
}


motif_coverage = function(){
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_motif_coverage_skip300.csv' 
  df = read_csv(infile)
  df = transform(df, unique_motif = as.factor(unique_motif))
  print(df$source)
  # stop()
  print(sapply(df, class))
  df$source[df$source == 'ab'] <- 'Paratope'
  df$source[df$source == 'ag'] = 'Epitope'
  df$source = factor(df$source, levels = c("Paratope", 'Epitope')) 
  medpara = aggregate(df[df$source=='Paratope',4], by= list(df[df$source=='Paratope',]$pdbidx),median)
  medepi = aggregate(df[df$source=='Epitope',4], by= list(df[df$source=='Epitope',]$pdbidx),median)
  print(head(df))
  meds  = df %>% group_by(pdbidx, source, strpdbidx) %>% summarise(med = round(median(coverage),2))
  print(tail(medpara))
  print(tail(medepi))
  print(meds)
  ggplot(df, aes(x=pdbidx, y=coverage)) + 
    # geom_point()+
    # geom_point(mapping = aes(x=pdbidx, y=coverage)) + 
    geom_boxplot(mapping = aes(x= reorder(strpdbidx, pdbidx), y=coverage*100, fill= source), lwd=0, outlier.shape = NA, flatten=0) +
    facet_wrap(~source, scales = 'free') +
    scale_fill_manual(values = c(abcolor3, agcolor), name = "Motif source", labels = c('Paratope', 'Epitope')) + 
    theme(axis.text.x = element_text(angle = 90)) +
    geom_text(data = meds, mapping = aes(x=strpdbidx, y=med*100+5, label = med), angle=90) +
    labs(x = 'Number of structures', y = 'Interaction motif coverage')
    # coord_flip()
    # geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE)
  outpdf(infile, 'motif_coverage', height = 10, width = 15)
}

branch_subset = function(){
  infile = 'abdb_outfiles_2019/paratope_epitope_internet_edges.csv'
  df = read_csv(infile)
  dfsub = df[df$source == 'XXXX1X', ]
  print(df)
  print(dfsub)
  print(sum(dfsub$weight))
}

# run stuff
# dl_ld_binary()
# chao_diversity()
# motif_coverage()
branch_subset()