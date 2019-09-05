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


abresnum_distribution_5a  <- function(infile){
# infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments.csv'
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
outname = sprintf('%s/%s_abresnum_distribution_sup3.pdf',outfigdir,inname)
ggsave(outname, height = 10, width = 8)
system(sprintf('open %s', outname))
###
}



ab.ag.venn = function(infile_ab, infile_ag){
  # infile_ab = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx.csv' 
  # infile_ag = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx.csv' 
  # infile_ag = 'abdb_outfiles/respairs_epitope_segment_notationx_top-1_gap_patterns.csv' 
  dfab = read_csv(infile_ab)
  dfag = read_csv(infile_ag)
  motif_ab = dfab$gap_patterns
  # print(motif_ab)
  motif_ag = dfag$gap_patterns
  x = list(motif_ab,motif_ag)
  partitions = get.venn.partitions(x)
  part_list = partitions$..values..
  percents = c(length(part_list$'1')/(length(part_list$'1') + length(part_list$'3')), length(part_list$'1')/(length(part_list$'1') + length(part_list$'2')))
  percents = round(percents,2)
  overlapdf = dfab[dfab$gap_patterns %in% part_list$'1', ]
  continuous_count = count(overlapdf, gapstrstatus)
  continuous_count$percent = round(continuous_count$n/(sum(continuous_count$n)),2)
  print(percents)
  print(continuous_count)
  print(infile_ab)
  inname = strsplit(tail(strsplit(infile_ab, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_venn_sup3.tiff', outfigdir,inname)
  # print(outname)
  names(x) = c('motif_ab', 'motif_ag')
  venn.plot <- venn.diagram(x=x, filename=outname,
                            output = FALSE ,
                            imagetype="tiff" ,
                            # height = 800 , 
                            # width = 900 , 
                            resolution = 300,
                            compression = "lzw",
                            lwd = 2,
                            lty = 'blank',
                            fill = c(abcolor3, agcolor),
                            alpha=c(0.7,0.7),
                            rotation.degree = 180,
                            cex = 5, 
                            cat.fontface=3,
                            cat.cex = 3,
                            # cat.dist = c(-0.05,-0.15),
                            cat.just = list(c(1.4, 0), c(-1.4, -3.4)),
                            category.names=c("Paratope", "Epitope")
                            )
  system(sprintf('open %s', outname))
  
}


# run stuff
# abresnum_distribution_5a('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments.csv')
# abresnum_distribution_5a('abdb_outfiles_2019/respairs_absort_cutoff4_abresnumi_segments.csv')
# abresnum_distribution_5a('abdb_outfiles_2019/respairs_absort_cutoff6_abresnumi_segments.csv')
ab.ag.venn(
  infile_ab = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx.csv', 
  infile_ag = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx.csv' 
)

ab.ag.venn(
  infile_ab = 'abdb_outfiles_2019/respairs_absort_cutoff6_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx.csv', 
  infile_ag = 'abdb_outfiles_2019/respairs_absort_cutoff6_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx.csv' 
)


ab.ag.venn(
  infile_ab = 'abdb_outfiles_2019/respairs_absort_cutoff4_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx.csv', 
  infile_ag = 'abdb_outfiles_2019/respairs_absort_cutoff4_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx.csv' 
)
