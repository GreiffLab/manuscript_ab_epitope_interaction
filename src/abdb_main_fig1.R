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
library('fossil')
library('gtools')

theme_set(theme.akbar())
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


abresnum.distribution  <- function(){
infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments.csv'
print(getwd())
df = read.csv(infile, stringsAsFactors = FALSE)
dt = setDT(df)
recounts_segment = dt[, .(.N), by = .(segment)]
total_res = sum(recounts_segment$N)
print(unique(df$segment))
respos = c(15, 33,57,78,98,38,108)
hrespos = list('HFR1'=15, 'CDR-H1'=33, 'CDR-H2'=57, 'HFR3'=78, 'CDR-H3' = 98,'HFR2'= 38, 'HFR4' = 108)
hchain_segment = recounts_segment[grep('H', segment)]
hchain_segment$fraction = hchain_segment$N/sum(hchain_segment$N)
respos = sapply(hchain_segment$segment, function(x) hrespos[x])
respos = unlist(respos)
# hchain_segment$respos = sapply(hchain_segment$segment, function(x) hrespos[x])
hchain_segment$respos = respos
lchain_segment = recounts_segment[grep('L', segment)]
lchain_segment$fraction = lchain_segment$N/sum(lchain_segment$N)
lrespos = list('LFR1'=11, 'CDR-L1'=28, 'CDR-L2'=54, 'LFR3'=68, 'CDR-L3' = 93,'LFR2'= 42, 'LFR4' = 105)
respos = unlist(sapply(lchain_segment$segment, function(x) lrespos[x]))
lchain_segment$respos = respos 
df[df=='H'] = 'Heavy'
df[df=='L'] = 'Light'
print(summary(df))
print(lchain_segment)
print(hchain_segment)
print(total_res)
xlab = sprintf('Total number of interacting residues: %s', total_res)
ggplot(data = df, aes(abchain, abresnum, fill=segment)) +
  geom_violin(fill = alpha('grey', 0.3), color=NA) +
  # geom_jitter(fill = alpha('red',0.3), color = alpha('red', 0.2)) +
  geom_jitter(shape=21, size=5, stroke=NA) +
  # geom_text(data=meds, mapping=aes(x=15, y=10, label = med_labels)) + 
  geom_text(data=hchain_segment, mapping=aes(x=1,y=respos, label= paste0(segment,': ', round(fraction*100,1), '%')), color='black') +
  geom_text(data=lchain_segment, mapping=aes(x=2,y=respos, label= paste0(segment,': ', round(fraction*100,1), '%')), color='black') +
  # scale_color_manual(values = my_spectral) +
  scale_fill_manual(values = cdrfr_colors) +
  labs(x=xlab, y='Residue number (Martin numbering scheme)', fill='Region') +
  theme(legend.text = element_text(size=8)) + 
  theme(legend.title = element_text(vjust = 1))
outname = 'abdb_figures/abresnum_distribution.pdf'
inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
outname = sprintf('%s/%s_abresnum_distribution.pdf',outfigdir,inname)
ggsave(outname, height = 10)
system(sprintf('open %s', outname))
###
}


length.distribution <- function(infile){
# infile = 'abdb_outfiles/respairs_paratope_segment.csv'
# infile = 'abdb_outfiles_2019/respairs_absort_cutoff6_abresnumi_segments_abshift_abshiftl_paratope_segment.csv'
inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1),'\\.')[[1]][1]
outfile = sprintf('%s/%s_length_distribution.pdf', outfigdir, inname)
print(outfile)
df = read.csv(infile, stringsAsFactors = FALSE)
if (grepl('epitope', infile)){
  df$plen = df$elen
} 
df[df == 'H'] <- 'Heavy chain'
df[df == 'L'] = 'Light chain'
dt = setDT(df)
print(dt)
meds = dt[, .(median(plen)),by = .(segment)]
maxs = dt[, .(max(plen)),by = .(segment)]
hmeds = meds[grep('H', segment)]
lmeds = meds[grep('L', segment)]
hmaxs = maxs[grep('H', segment)]
lmaxs = maxs[grep('L', segment)]
colnames(hmeds) = c('segment', 'mean')
colnames(hmaxs) = c('segment', 'max')
colnames(lmeds) = c('segment', 'mean')
colnames(lmaxs) = c('segment', 'max')
print(summary(df))
# beeswarm(segment ~ plen, data=df )
  ggplot(data=df) +
  geom_point(aes(x=segment, y=plen, fill=segment), size = 2, shape = 21, stroke = NA, position = position_jitterdodge(jitter.height = 0.2, jitter.width = 0.5, dodge.width = 0.5)) +
  geom_boxplot(mapping = aes(x=segment, y=plen, fill=segment), outlier.shape = NA, lwd=0, fatten=0) +
  geom_text(data=hmeds, size =4, mapping = aes(x=segment, y= mean+0.3, label = paste0('median: ', round(mean,1)))) +
  geom_text(data=hmaxs, size =4, mapping = aes(x=segment, y= max+1, label = paste0('max: ', round(max,1)))) +
  geom_text(data=lmeds, size=4, mapping = aes(x=segment, y=mean+0.3, label = paste0('median: ', round(mean,1)))) +
  geom_text(data=lmaxs, size =4, mapping = aes(x=segment, y= max+1, label = paste0('max: ', round(max,1)))) +
  labs(x='Region', y = '# of interacting residues', fill='Chain') +
  # facet_wrap(~ abchain, scale = 'free') + 
  scale_color_manual(values = my_spectral) + 
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title = element_text(size=25)) + 
  theme(legend.position = 'none') + 
  scale_fill_manual(values = cdrfr_colors) 
# outname = 'abdb_figures/paratope_length_distribution.pdf'
ggsave(outfile, width=15)
system(sprintf('open %s', outfile ))
###
}

gap.distribution <- function(infile){
  # infile = 'abdb_outfiles/respairs_paratope_segment_gaponly.csv'
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_gaponly.csv'
  df = read.csv(infile, stringsAsFactors = FALSE)
  df = subset(df, (gap < 1000)&(gap>0))
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outfile = sprintf('%s/%s_gap_distribution.pdf', outfigdir, inname)
  print(head(df))
  if (grepl('paratope', infile)){
    maxadj = 1
    medadj = 0.4
  } else if (grepl('epitope', infile)){
    maxadj = 20
    medadj = 7
  } 
  dt = setDT(df)
  print(dt)
  print(dt[, segment=='CDR-H1' & gap > 0, gap])
  meds = dt[, .(median(gap)),by = .(segment)]
  maxs = dt[, .(max(gap)),by = .(segment)]
  hmeds = meds[grep('H', segment)]
  lmeds = meds[grep('L', segment)]
  hmaxs = maxs[grep('H', segment)]
  lmaxs = maxs[grep('L', segment)]
  colnames(hmeds) = c('segment', 'mean')
  colnames(hmaxs) = c('segment', 'max')
  colnames(lmeds) = c('segment', 'mean')
  colnames(lmaxs) = c('segment', 'max')
  xorder =  c('CDR-H1', 'CDR-H2', 'CDR-H3', 'HFR1', 'HFR2', 'HFR3', 'HFR4', 'CDR-L1', 'CDR-L2', 'CDR-L3', 'LFR1', 'LFR2', 'LFR3')
  # stop()
  df[df == 'H'] <- 'Heavy chain'
  df[df == 'L'] <- 'Light chain'
  ggplot(data=df) +
  geom_point(aes(x=segment, y=gap, fill=segment), size = 2, shape = 21, stroke = NA, position = position_jitterdodge(jitter.height = 0.2, jitter.width = 0.5, dodge.width = 0.5)) +
  geom_boxplot(mapping = aes(x=segment, y=gap, fill=segment), outlier.shape = NA, lwd=0, fatten=0) +
  geom_text(data=hmeds, size =4, mapping = aes(x=segment, y= mean+medadj, label = paste0('median: ', round(mean,1)))) +
  geom_text(data=hmaxs, size =4, mapping = aes(x=segment, y= max+maxadj, label = paste0('max: ', round(max,1)))) +
  geom_text(data=lmeds, size=4, mapping = aes(x=segment, y=mean+medadj, label = paste0('median: ', round(mean,1)))) +
  geom_text(data=lmaxs, size =4, mapping = aes(x=segment, y= max+maxadj, label = paste0('max: ', round(max,1)))) +
  # facet_wrap(~ abchain, scales = 'free_x') +
  scale_fill_manual(values = cdrfr_colors) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title = element_text(size=25)) + 
  labs(x='Region', y= 'Gap length (# of non-interacting residues)') + 
  theme(legend.position = 'none')
  # outname = 'abdb_figures/paratope_gap_distribution.pdf'
  ggsave(outfile, width=15)
  system(sprintf('open %s', outfile))
  print(summary(df))  
}



# run stuff
# abresnum.distribution()
leninfile_paratope = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment.csv'
leninfile_epitope = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment.csv'
length.distribution(leninfile_paratope)
length.distribution(leninfile_epitope)
gapinfile_paratope = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_gaponly.csv'
gapinfile_epitope = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_gaponly.csv'
gap.distribution(gapinfile_paratope)
gap.distribution(gapinfile_epitope)

