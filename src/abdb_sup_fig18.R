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
library('wesanderson')
names(wes_palettes)
theme_set(themeakbar())
my_spectral <- colorRampPalette(brewer.pal(8,'Spectral'))(14)
my_darjee = colorRampPalette(wes_palette('Darjeeling1'))(8)
print(my_darjee)
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

my_ccols = c('salmon',
             'deepskyblue1', 'deepskyblue2', 'deepskyblue3', 'deepskyblue4',
             'orange1', 'orange2', 'orange3', 'orange4')


my_ccols2 = c('salmon',
             abcolor3, abcolor3, abcolor3, abcolor3,
             ppi_mcolor, ppi_mcolor, ppi_mcolor, ppi_mcolor)

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

gap_status = function(x){
  x[x > 0] = 1
  x[x == 1] = 'discontinuous'
  x[x == 0] = 'continuous'
  return = x
}

pattern.bar.pie <- function(infile){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv'
  content <- read_csv(infile)
  print(content)
  content1 = content %>% 
    mutate(gap_status1 = gap_status(max_gap1)) %>%
    count(., gap_pattern1, gap_status1, intertype1)
  print(content1)
  content2 = content %>% 
    mutate(gap_status2 = gap_status(max_gap2)) %>%
    count(., gap_pattern2, gap_status2, intertype1)
  print(content2)
  g <- ggplot(data=content1)
  g <- g + geom_bar(stat = 'identity', mapping = aes(x='', y= n, fill = gap_status1, color = gap_status1))
  g <- g + facet_wrap(~ intertype1, ncol = 3)
  g <- g + coord_polar('y', start=0)
  g <- g + theme(axis.text.y = element_blank(), 
                 axis.title.x = element_text(size=20),
                 strip.text = element_text(size = 20)) +
    labs(x=' ', y='Number of continuous and discontinous motifs') +
    scale_fill_manual(values=c('red', 'black')) +
    scale_color_manual(values=c('red', 'black')) +
    labs (fill= 'Gap Status', color= 'Gap Status') + 
    theme(legend.position = 'none')
  outpng(infile, 'motif', height = 5)
  ### motif partner
  g <- ggplot(data=content2)
  g <- g + geom_bar(stat = 'identity', mapping = aes(x='', y= n, fill = gap_status2, color = gap_status2))
  g <- g + facet_wrap(~ intertype1, ncol = 3)
  g <- g + coord_polar('y', start=0)
  g <- g + theme(axis.text.y = element_blank(),
                 axis.title.x = element_text(size=20),
                 strip.text = element_text(size = 20)) +
    labs(x=' ', y='Number of continuous and discontinous motifs') +
    scale_fill_manual(values=c('red', 'black')) +
    scale_color_manual(values=c('red', 'black')) +
    labs (fill= 'Gap Status', color= 'Gap Status') + 
    theme(legend.position = 'none')
  outpng(infile, 'motif_partner', height = 5)
    
}


pattern.bar.radial <- function(infile){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv'
  content <- read_csv(infile)
  print(content)
  content1 = content %>% 
    mutate(gap_status1 = gap_status(max_gap1)) %>%
    count(., gap_pattern1, gap_status1, intertype1)
  print(content1)
  content2 = content %>% 
    mutate(gap_status2 = gap_status(max_gap2)) %>%
    count(., gap_pattern2, gap_status2, intertype1)
  print(content2)
  ggplot(data=content1) +
    geom_bar(stat = 'identity', mapping = aes(x=gap_pattern1, y= n, fill = gap_status1, color = gap_status1)) + 
    facet_wrap(~ intertype1, ncol = 3) + 
    coord_polar('y', start=0) +
    labs(x='Motif', y='Number of shared motifs') +
    theme(axis.title = element_text(size = 20), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(),
          strip.text = element_text(size = 20)) +
    scale_fill_manual(values=c('red', 'black')) +
    scale_color_manual(values=c('red', 'black')) +
    labs (fill= 'Gap Status', color= 'Gap Status') + 
    theme(legend.position = 'none') + 
    scale_y_continuous(breaks = seq(0,1600, by = 500))
  outpng(infile, 'radial_motif', height = 5)
  # stop(000)
  ### motif partner
  g <- ggplot(data=content2)
  g <- g + geom_bar(stat = 'identity', mapping = aes(x=gap_pattern2, y= n, fill = gap_status2, color = gap_status2))
  g <- g + facet_wrap(~ intertype1, ncol = 3) + 
    coord_polar('y', start=0) +
    labs(x='Motif', y='Number of shared motifs') +
    labs(x='Motif', y='Number of shared motifs') +
    theme(axis.title = element_text(size = 20), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(),
          strip.text = element_text(size = 20)) +
    scale_fill_manual(values=c('red', 'black')) +
    scale_color_manual(values=c('red', 'black')) +
    labs (fill= 'Gap Status', color= 'Gap Status') + 
    theme(legend.position = 'none') +
    scale_y_continuous(breaks = seq(0,1600, by = 500))
  outpng(infile, 'radial_motif_partner', height = 5)
    
}



pattern.bar.radial2 <- function(infile){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_count.csv'
  outname = strsplit(infile, split ='/')[[1]][2]
  outname = paste0('abdb_figures_2019/',substr(outname, 1,nchar(outname)-4), '_radial.pdf') 
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
    theme(legend.position = 'none')
  ggsave(outname, width=8, height=15)
  system(sprintf('open %s', outname))
}


shared_structures = function(){
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv'
  df = read_csv(infile)
  countdf  = count(df, pdbid1, gap_pattern1)
  countdf = count(countdf, gap_pattern1)
  countdf2  = count(df, pdbid1, gap_pattern2)
  countdf2 = count(countdf2, gap_pattern2)
  print(countdf2)
  outdf = data_frame()
  for (i in seq(2,20)){
    print(i)
    sdf = countdf[countdf$n>=i,]
    # print(dim(sdf))
    newrow = data_frame(nmotif = dim(sdf)[1], cutoff = i, source = 'PPI motif', motif_percent = dim(sdf)[1]/dim(countdf)[1])
    outdf = rbind(outdf, newrow)
    sdf2 = countdf2[countdf2$n>i,]
    # print(dim(sdf))
    newrow2 = data_frame(nmotif = dim(sdf2)[1], cutoff = i, source = 'PPI partner', motif_percent = dim(sdf2)[1]/dim(countdf2)[1])
    outdf = rbind(outdf, newrow2)
    # stop()
  }
  print(outdf)
  outdf$source = factor(outdf$source, level = c('PPI motif', 'PPI partner'))
  ggplot(data=outdf) +
    geom_bar(mapping = aes(x=cutoff, y=nmotif, fill=source), stat = 'identity') + 
    facet_wrap(~source, scales = 'free')+ 
    geom_text(data=outdf, mapping = aes(x=cutoff, y=nmotif, label = sprintf('%s\n%s%%', nmotif, round(motif_percent,2)*100)), nudge_y = 40) + 
    # labs(x= 'Cut-offs for number of structures \n(2-20 structures)', y = 'Motifs found in a cut-off') + 
    labs(x= '# of PPI complexes\n(2-20 structures)', y = '# of motifs found across a given\nPPI complexes') + 
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 25),
          legend.title = element_blank(), 
          # legend.text = element_text(size=25),
          legend.position = 0,
          strip.text = element_text(size = 25),) +
    scale_fill_manual(values = c(ppi_mcolor, ppi_pcolor))
  outpdf(infile, 'ppi_shared_structures', width = 20)
}




# run stuff
# chao_diversity()
# motif_coverage()
# pattern.bar.pie('abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv')
# pattern.bar.radial('abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv')
shared_structures()

