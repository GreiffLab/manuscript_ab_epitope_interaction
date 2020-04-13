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
print(my_spectral)
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



shm_resnum_distribution = function(cutoff){
  infile = 'abdb_outfiles_2019/homo_mus_shm_residues.csv'
  # infile = 'abdb_outfiles_2019/shm_residues.csv'
  df = read_csv(infile)
  print(df)
  print(dim(df))
  infile2 = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_adjusted_resnum.csv'
  df2 = read_csv(infile2)
  df = df[df$identity > cutoff,]
  print(df)
  df[df == 'H'] = 'Heavy'
  df[df == 'L'] = 'Light'
  df2[df2 == 'H'] = 'Heavy'
  df2[df2 == 'L'] = 'Light'
  df$resnum_adj = as.character(as.numeric(df$resnum)+1)
  countdf2 = drop_na(count(df, resnum_adj, chain))
  countdf = count(df, chain_gene, pdbid)
  chaincountdf = count(df, chain, pdbid) %>% count(chain)
  chaincount_label  = sprintf('Total number of structures: %s', chaincountdf$n)
  rescountdf = count(df, resnum_adj, chain, sort = TRUE)
  rescountdf = rescountdf[order(rescountdf$resnum_adj),]
  rescountdf2 = count(df2, adj_resnum, abchain, sort =TRUE)
  rescountdf2$adj_resnumstr = as.character(rescountdf2$adj_resnum)
  rescountdf2 = rescountdf2[order(rescountdf2$adj_resnumstr),]
  print(rescountdf)
  print(rescountdf2)
  # write.csv(rescountdf, 'abdb_outfiles_2019/shm_paratope_adjusted_rescountdf.csv', row.names = FALSE)
  # write.csv(rescountdf2, 'abdb_outfiles_2019/shm_paratope_adjusted_rescountdf2.csv', row.names = FALSE)
  #continue in python. need to do some gymnastic to account for non existant residues in paratope
  # stop()
  # rescountdf2$n_shm = rescountdf$n[1:dim(rescountdf2)[1]]
  heavy_rescountdf2 = rescountdf2[rescountdf2$abchain=='Heavy',]
  heavy_rescountdf = rescountdf[rescountdf$chain=='Heavy',]
  heavy_rescountdf = heavy_rescountdf[heavy_rescountdf$resnum_adj %in%heavy_rescountdf2$adj_resnum,]
  heavy_rescountdf2 = heavy_rescountdf2[heavy_rescountdf2$adj_resnum %in%heavy_rescountdf$resnum_adj,]
  print(heavy_rescountdf)
  print(heavy_rescountdf2)
  light_rescountdf2 = rescountdf2[rescountdf2$abchain=='Light',]
  light_rescountdf = rescountdf[rescountdf$chain=='Light',]
  light_rescountdf = light_rescountdf[light_rescountdf$resnum_adj %in%light_rescountdf2$adj_resnum,]
  light_rescountdf2 = light_rescountdf2[light_rescountdf2$adj_resnum %in%light_rescountdf$resnum_adj,]
  pcor_heavy = cor.test(heavy_rescountdf$n, heavy_rescountdf2$n)
  pcor_light = cor.test(light_rescountdf$n, light_rescountdf2$n)
  scor_heavy = cor.test(heavy_rescountdf$n, heavy_rescountdf2$n, method = 'spearman')
  scor_light = cor.test(light_rescountdf$n, light_rescountdf2$n, method = 'spearman')
  pcordf = as.data.frame(list('chain' = c('Heavy', 'Light'), 'cor' = c(pcor_heavy$estimate, pcor_light$estimate),
                         'pval' = c(pcor_heavy$p.value, pcor_light$p.value)))
  scordf = as.data.frame(list('chain' = c('Heavy', 'Light'), 'cor' = c(scor_heavy$estimate, scor_light$estimate),
                         'pval' = c(scor_heavy$p.value, scor_light$p.value)))
  # pcordf$pval[pcordf$pval < 10e-10] = '<10e-10'
  # scordf$pval[scordf$pval < 10e-10] = '<10e-10'
  print(pcordf)
  pcordf$pval = signif(pcordf$pval, digits = 2)
  scordf$pval = signif(scordf$pval, digits = 2)
  # print(pcordf)
  # print(scordf)
  # stop()
  pcorlabel = sprintf('Correlation to paratope (Pearson): %s, p-value: %s', round(pcordf$cor,2), pcordf$pval)
  scorlabel = sprintf('Correlation to paratope (Spearman): %s, p-value %s', round(scordf$cor,2), scordf$pval)
  ggplot(data = countdf2) + 
    geom_bar(stat = 'identity', mapping = aes(x=reorder(resnum_adj, as.numeric(resnum_adj)), y=n), color=NA, fill = alpha(my_spectral[12], 0.8)) +
    facet_wrap(~ chain) +
    theme(axis.text.x = element_text(angle = 90)) + 
    labs(x=' ', y=' ') +
    xlim(as.character(seq(1,131))) +
    geom_text(data=pcordf, mapping = aes(x=5, y = 300, label=pcorlabel), size = 7, hjust = 0) + 
    geom_text(data=pcordf, mapping = aes(x=5, y = 300-20, label=scorlabel), size=7, hjust = 0) + 
    geom_text(data=chaincountdf, mapping = aes(x=5, y = 300-40, label=chaincount_label), size=7, hjust = 0) + 
    theme(legend.position = 0,
          strip.text = element_text(size = 30),
          axis.text.y = element_text(size = 25)) +
    ylim(c(0,300))
  outpdf(infile, sprintf('resnum_distribution_id%s', cutoff*100), width = 35)
}


paratope_resnumadj_distribution = function(){
  infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_adjusted_resnum.csv'
  df = read_csv(infile)
  print(df)
  df[df=='H'] = 'Heavy'
  df[df=='L'] = 'Light'
  df$adj_resnum = as.character(df$adj_resnum)
  countdf2 = drop_na(count(df, adj_resnum, abchain))
  ggplot(data = countdf2) + 
    geom_bar(stat = 'identity', mapping = aes(x=reorder(adj_resnum, as.numeric(adj_resnum)), y=n), color=NA, fill = abcolor3) +
    facet_wrap(~abchain) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    labs(x='Residue position', y=' ') +
    theme(legend.position = 0,
          strip.text = element_text(size = 30),
          axis.text.y = element_text(size = 25)) + 
    xlim(as.character(seq(1,131)))
  outpdf(infile,'adj_resnum_distribution', width = 35)
}

shm_paratope_correlation = function(){
  infile = 'abdb_outfiles_2019/shm_residues.csv'
  infile2 = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_adjusted_resnum.csv'
  df = read_csv(infile)
  df2 = read_csv(infile2)
  print(df)
  print(df2)

                 
}


shm_paratope_resnum_cor = function(cutoff){
  infile = 'abdb_outfiles_2019/shm_residues.csv'
  df = read_csv(infile)
  infile2 = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_adjusted_resnum.csv'
  df2 = read_csv(infile2)
  df = df[df$identity > cutoff,]
  print(df)
  df[df == 'H'] = 'Heavy'
  df[df == 'L'] = 'Light'
  df2[df2 == 'H'] = 'Heavy'
  df2[df2 == 'L'] = 'Light'
  print(df2)
  # stop()
  df$resnum_adj = as.character(as.numeric(df$resnum)+1)
  countdf2 = drop_na(count(df, resnum_adj, chain))
  # print(countdf)
  # print(df[df$res == '-',])
  countdf = count(df, chain_gene, pdbid)
  rescountdf = count(df, resnum_adj, chain, sort = TRUE)
  rescountdf = rescountdf[order(rescountdf$resnum_adj),]
  rescountdf2 = count(df2, adj_resnum, abchain, sort =TRUE)
  rescountdf2$adj_resnumstr = as.character(rescountdf2$adj_resnum)
  rescountdf2 = rescountdf2[order(rescountdf2$adj_resnumstr),]
  print(rescountdf)
  print(rescountdf2)
  # write.csv(rescountdf, 'abdb_outfiles_2019/shm_paratope_adjusted_rescountdf.csv', row.names = FALSE)
  # write.csv(rescountdf2, 'abdb_outfiles_2019/shm_paratope_adjusted_rescountdf2.csv', row.names = FALSE)
  #continue in python. need to do some gymnastic to account for non existant residues in paratope
  # stop()
  # rescountdf2$n_shm = rescountdf$n[1:dim(rescountdf2)[1]]
  heavy_rescountdf2 = rescountdf2[rescountdf2$abchain=='Heavy',]
  heavy_rescountdf = rescountdf[rescountdf$chain=='Heavy',]
  heavy_rescountdf = heavy_rescountdf[heavy_rescountdf$resnum_adj %in%heavy_rescountdf2$adj_resnum,]
  heavy_rescountdf2 = heavy_rescountdf2[heavy_rescountdf2$adj_resnum %in%heavy_rescountdf$resnum_adj,]
  print(heavy_rescountdf)
  print(heavy_rescountdf2)
  light_rescountdf2 = rescountdf2[rescountdf2$abchain=='Light',]
  light_rescountdf = rescountdf[rescountdf$chain=='Light',]
  light_rescountdf = light_rescountdf[light_rescountdf$resnum_adj %in%light_rescountdf2$adj_resnum,]
  light_rescountdf2 = light_rescountdf2[light_rescountdf2$adj_resnum %in%light_rescountdf$resnum_adj,]
  pcor_heavy = round(cor(heavy_rescountdf$n, heavy_rescountdf2$n),2)
  pcor_light = round(cor(light_rescountdf$n, light_rescountdf2$n),2)
  scor_heavy = round(cor(heavy_rescountdf$n, heavy_rescountdf2$n, method = 'spearman'),2)
  scor_light = round(cor(light_rescountdf$n, light_rescountdf2$n, method = 'spearman'),2)
  pcordf = as.data.frame(list('chain' = c('Heavy', 'Light'), 'cor' = c(pcor_heavy, pcor_light)))
  scordf = as.data.frame(list('chain' = c('Heavy', 'Light'), 'cor' = c(scor_heavy, scor_light)))
  print(pcor_heavy)
  print(pcordf)
  heavy_rescountdf2$n_shm = heavy_rescountdf$n
  light_rescountdf2$n_shm = light_rescountdf$n
  print(heavy_rescountdf2)
  mergeddf = rbind(heavy_rescountdf2,light_rescountdf2)
  print(mergeddf)
  pcorlabel = sprintf('Correlation to paratope (Pearson): %s', pcordf$cor)
  scorlabel = sprintf('Correlation to paratope (Spearman): %s', scordf$cor)
  ggplot(data = mergeddf) + 
    # geom_bar(stat = 'identity', mapping = aes(x=reorder(resnum_adj, as.numeric(resnum_adj)), y=n), color=NA, fill = alpha(my_spectral[12], 0.8)) +
    geom_point(mapping = aes(x=n, y = n_shm), size=7, color = alpha(abcolor3,0.5)) +
    facet_wrap(~ chain) +
    theme(axis.text.x = element_text(angle = 90)) + 
    labs(x='Paratope', y='SHM') +
    # xlim(as.character(seq(1,131))) +
    geom_text(data=pcordf, mapping = aes(x=250, y = max(countdf2$n), label=pcorlabel), size = 7) + 
    geom_text(data=pcordf, mapping = aes(x=250, y = max(countdf2$n)-15, label=scorlabel), size=7) + 
    theme(legend.position = 0)
  outpdf(infile, sprintf('resnum_correlation_id%s', cutoff*100), width = 35)
  write.csv(rescountdf2, 'abdb_outfiles_2019/shm_paratope_adjustedres_count.csv', row.names = FALSE)
}


# run stuff
shm_resnum_distribution(0.7)
shm_resnum_distribution(0.8)
shm_resnum_distribution(0.9)
# paratope_resnumadj_distribution()
# shm_paratope_resnum_cor(0.9)
# shm_paratope_resnum_cor(0.8)
# shm_paratope_resnum_cor(0.7)

# heavy light. align index, highlight cdrs. number complexes per cutoff. drop the inferred. number of complexes with mutations.