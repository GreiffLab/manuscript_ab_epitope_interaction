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


ppi_abdb_venn = function(){
  infile_ab = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx.csv' 
  infile_ag = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx.csv' 
  infile_ppi = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300.csv'
  infilemerged = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil_pc.csv'
  mergedf  = read_csv(infilemerged)
  print(mergedf)
  # stop()
  abdf = read_csv(infile_ab)
  abdf = abdf[!is.na(abdf$paratope),]
  abmotif = abdf$gap_patterns
  agdf = read_csv(infile_ag)
  agdf = agdf[!is.na(agdf$epitope),]
  agmotif = agdf$gap_patterns
  ppidf = read_csv(infile_ppi)
  print(ppidf)
  print(abdf)
  print(agdf)
  # stop()
  nlabels = sprintf('Total para: %s \nTotal epi: %s \nTotal PPI: %s', dim(mergedf)[1], dim(mergedf)[1], dim(ppidf)[1])
  print(nlabels)
  punique = length(unique(mergedf$ab_motif))
  eunique = length(unique(mergedf$ag_motif))
  ppiunique = length(unique(ppidf$gap_pattern))
  uniquelabel = sprintf('Unique para: %s \nUnique epi: %s \nUnique PPI: %s', punique, eunique, ppiunique)
  print(uniquelabel)
  # stop()
  # ppidf = subset(ppidf, max_gap < 8)
  ppimotif = ppidf$gap_pattern
  x = list(abmotif, agmotif, ppimotif)
  outname  = 'abdb_figures_2019/paratope_epitope_ppi_venn.tiff'
  venn.plot <- venn.diagram(x=x, filename=outname,
                            output = FALSE ,
                            imagetype="tiff" ,
                            # height = 800 , 
                            # width = 900 , 
                            resolution = 300,
                            compression = "lzw",
                            lwd = 2,
                            lty = 'blank',
                            fill = c(abcolor3, agcolor, 'orange'),
                            alpha=c(0.5,0.5, 0.7),
                            cex = 5, 
                            cat.fontface=1,
                            cat.cex = 3,
                            cat.dist = c(-0.004,-0.001, 0.0001),
                            category.names=c("Paratope", "Epitope", "PPI"),
                            main = nlabels,
                            main.cex = 2,
                            main.pos = c(0.15,0.15),
                            sub = uniquelabel,
                            sub.cex = 2,
                            sub.pos = c(0.85,0.15)
  ) 
  system(sprintf('open %s', outname))
}

ppi_motif_len = function(){
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300.csv'
  df = read_csv(infile)
  print(df)
  # stop()
  # df = df[seq(100),]
  # df = subset(df, (max_gap < 8))
  medians = aggregate(df, by = list(df$intertype), median)
  medians = df %>% group_by(intertype) %>% summarise(max_len = max(motif_len), med_len = median(motif_len))
  print(medians)
  ggplot(data= df, mapping = aes(x=intertype, y=motif_len, fill=intertype, color = intertype)) + 
    geom_violin(color=NA) + 
    # geom_jitter(color = alpha('orange', 0.01), size=5) + 
    scale_fill_manual(values = c(alpha('orange'), alpha('orange',0.5))) +
    scale_color_manual(values = c(alpha('orange'), alpha('orange',0.5))) +
    geom_text(data=medians, mapping = aes(x=intertype, y=med_len), label = sprintf('Median: %s', medians$med_len), color = 'Black') +
    geom_text(data=medians, mapping = aes(x=intertype, y=max_len), label = sprintf('Max: %s', medians$max_len), color = 'Black') +
    labs(x= 'Domain type', y = 'PPI interaction motif length') + 
    theme(legend.position = 'none', 
          axis.title = element_text(size = 25),
          axis.text = element_text(size = 25))
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_ppi_motif_len.pdf', outfigdir,inname)
  ggsave(outname, width = 6, height = 9)
  system(sprintf('open %s', outname))
}


ppi_max_gap = function(){
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300.csv'
  df = read_csv(infile)
  df = subset(df, max_gap < 8)
  # df = df[seq(100),]
  maxes = df %>% group_by(intertype) %>% summarise(ma_gap = max(max_gap), med_gap = median(max_gap))
  print(maxes)
  print(df$max_gap)
  # stop()
  ggplot(data= df, mapping = aes(x=intertype, y=max_gap, fill=intertype, color = intertype)) + 
    geom_violin() +
    # geom_jitter(color = alpha('orange', 0.01), size=5) + 
    scale_fill_manual(values = c(alpha('orange'), alpha('orange',0.5))) + 
    scale_color_manual(values = c(alpha('orange'), alpha('orange',0.5))) + 
    geom_text(data=maxes, mapping = aes(x=intertype, y=med_gap), label = sprintf('Median: %s', maxes$med_gap), color= 'Black') +
    geom_text(data=maxes, mapping = aes(x=intertype, y=ma_gap), label = sprintf('Median: %s', maxes$ma_gap), color= 'Black') +
    labs(x= 'Domain type', y = 'PPI maximum gap') + 
    theme(legend.position = 'none', 
          axis.title = element_text(size = 25),
          axis.text = element_text(size = 25))
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_ppi_max_gap.pdf', outfigdir,inname)
  ggsave(outname, width = 6, height = 9)
  system(sprintf('open %s', outname))
}

outpdf = function(infile, tag, width=8, height=8){
  # assumes infile has extensins (.csv etc)
  # uses the first bit of the file as name and add the tag
  # uses ggsave, opens pdf plot after saving
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_%s.pdf', outfigdir,inname,tag)
  ggsave(outname, width=width, height = height)
  print(outname)
  system(sprintf('open %s', outname))
  
}

ppi_residue_distribution = function(){
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_three_letter_residue.csv'
  df = read_csv(infile)
  print(df)
  ylab = sprintf('Percentage (%%) \n[Total number of interacting residues: %s]', dim(df)[1])
  print(ylab)
  # stop()
  getpdbid = function(input){
    print(str(input))
    inds = c(1:length(input))
    print(str(inds))
    print(substr(input[2], 1, 4))
    out  = sapply(inds, function(x) substr(input[x],1,4))
  }
  uniquepdb = df %>% mutate(pdbid = getpdbid(.$pdbchainpair)) %>% group_by(domain_type) %>% summarise(npdb = length(unique(pdbid)))
  print(uniquepdb)
  infile2 = 'abdb_outfiles_2019/threedid_no_iglike.csv'
  # no iglike numbers
  df2 = read_csv(infile2)
  uniquepdb2 = df2 %>% group_by(intertype) %>% summarise(npdb = length(unique(pdbid)))
  print(uniquepdb2)
  print('no iglike # of proteins')
  # stop()
  countdf = count(df, residue, domain_type) %>% group_by(domain_type) %>% mutate(percent = round(n/sum(n)*100,1))
  sums = countdf %>% group_by(domain_type) %>% summarise(sum = sum(percent))
  print(sums)
  ggplot(data=countdf) + 
    geom_bar(mapping = aes(x=residue, y = percent,fill= domain_type), stat='identity') + 
    facet_wrap(~domain_type) +
    scale_fill_manual(values = c(alpha('orange'), alpha('orange',0.5))) +
    geom_text(data= countdf, mapping = aes(x=residue, y=percent+0.5, label = percent)) + 
    geom_text(data= uniquepdb, mapping = aes(x=20.4, y=7.5, label = sprintf('# of proteins: %s', npdb)), hjust=0) + 
    # scale_y_continuous(labels= percent) + 
    theme(legend.position = 0, 
          axis.title = element_text(size = 30),
          axis.text = element_text(size = 20),
          strip.text = element_text(size = 20)) + 
    coord_flip() +
    labs(y=ylab, x = 'Residue name')
  outpdf(infile, 'box_lot', width = 11, height = 11)
}


ppi_residue_correlation = function(){
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_three_letter_residue.csv'
  df = read_csv(infile)
  print(df)
  countdf = count(df, residue, domain_type)
  countdf[countdf=='interdomain'] = 'PPI\ninter'
  countdf[countdf=='intradomain'] = 'PPI\nintra'
  print(countdf)
  infile2 = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments.csv'
  df2 = read_csv(infile2)
  countdf2 = count(df2, abres, abchain)
  countdf2[countdf2=='H'] = 'Paratope\nheavy'
  countdf2[countdf2=='L'] = 'Paratope\nlight'
  print(countdf2)
  setnames(countdf2, old = c('abres', 'abchain', 'n'), new=c('residue', 'domain_type', 'n'))
  countdf3 = count(df2, agres, abchain)
  countdf3[countdf3=='H'] = 'Epitope\nheavy'
  countdf3[countdf3=='L'] = 'Epitope\nlight'
  setnames(countdf3, old = c('agres', 'abchain', 'n'), new=c('residue', 'domain_type', 'n'))
  print(countdf3)
  mergedf = rbind(countdf, countdf2, countdf3)
  spreaddf = spread(mergedf, domain_type, n)
  cormat_data = spreaddf[, 2:dim(spreaddf)[2]]
  print(cormat_data)
  cormat = round(cor(cormat_data),2)
  cormat_spearman = round(cor(cormat_data, method = 'spearman'), 2)
  cormat[lower.tri(cormat)] = cormat_spearman[lower.tri(cormat_spearman)]
  print(cormat_spearman)
  print(cormat)
  cormat_melted = melt(cormat)
  print(cormat_melted)
  ggplot(data = cormat_melted, aes(Var2, Var1, fill = round(value,2)))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = abcolor3, high = my_spectral[3], mid = "white", 
                         midpoint = 0., limit = c(-0.05,1), space = "Lab", 
                         name="Upper: Spearman\nLower: Pearson") +
    theme_minimal()+ 
    theme(axis.text = element_text(angle = 0, vjust = 1,size = 20, hjust = 0.5),
          legend.position = 'bottom',
          legend.text = element_text(size = 9),
          legend) +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 10) +
    labs(x='', y = '')
    # theme(
    #   axis.title.x = element_blank(),
    #   axis.title.y = element_blank(),
    #   panel.grid.major = element_blank(),
    #   panel.border = element_blank(),
    #   panel.background = element_blank(),
    #   axis.ticks = element_blank(),
    #   legend.justification = c(1, 0),
    #   legend.position = c(0.6, 0.7),
    #   legend.direction = "horizontal")+
    # guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
    #                              title.position = "top", title.hjust = 0.5)) 
  outpdf(infile, 'cormat_plot', width = 11, height = 9)
}

ppi_motif_distribution_interdomain = function(){
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300.csv'
  df = read_csv(infile)
  print(df)
  df = df[df$intertype == 'interdomain',]
  print(df)
  countdf = count(df, gap_pattern, sort = TRUE)
  topndf = countdf[1:200,]
  print(topndf)
  ggplot(data=topndf) +
    geom_bar(mapping = aes(x=reorder(gap_pattern, -n), y=n), fill='orange', stat = 'identity') +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y='# of motifs', x = 'PPI motifs') +
    theme(axis.title = element_text(size = 60), legend.title = element_text(size = 50),
          legend.text = element_text(size=30)) +
    geom_text(mapping = aes(x= gap_pattern, y=n+80, label = n), angle=90 )
  outpdf(infile, 'motif_distribution', width = 32, height = 28)
  cumdf = topndf %>% mutate(csum = cumsum(n)) %>% mutate(csum_percent = percent(csum/sum(countdf$n)))
  print(cumdf)
  ggplot(data=cumdf) + 
    geom_bar(mapping = aes(x= reorder(gap_pattern, csum), y=csum), fill='orange', stat = 'identity') +
    theme(axis.text.x = element_blank()) +
    labs(y='Cummulative percentage (%)', x = '') + 
    theme(axis.title = element_text(size = 60), legend.title = element_text(size = 50),
          legend.text = element_text(size=30)) + 
    geom_text(mapping = aes(x= gap_pattern, y=csum+900, label = csum_percent), angle=90 )
  outpdf(infile, 'motif_csum', width = 32, height = 20)
  
  
  
}


ppi_motif_distribution_intradomain = function(){
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300.csv'
  df = read_csv(infile)
  print(df)
  df = df[df$intertype == 'intradomain',]
  print(df)
  countdf = count(df, gap_pattern, sort = TRUE)
  topndf = countdf[1:200,]
  print(topndf)
  ggplot(data=topndf) +
    geom_bar(mapping = aes(x=reorder(gap_pattern, -n), y=n), fill=alpha('orange',0.5), stat = 'identity') +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y='# of motifs', x = 'PPI motifs') +
    theme(axis.title = element_text(size = 60), legend.title = element_text(size = 50),
          legend.text = element_text(size=30)) +
    geom_text(mapping = aes(x= gap_pattern, y=n+8, label = n), angle=90 )
  outpdf(infile, 'motif_distribution_intra', width = 32, height = 28)
  cumdf = topndf %>% mutate(csum = cumsum(n)) %>% mutate(csum_percent = percent(csum/sum(countdf$n)))
  print(cumdf)
  ggplot(data=cumdf) + 
    geom_bar(mapping = aes(x= reorder(gap_pattern, csum), y=csum), fill= alpha('orange', 0.5), stat = 'identity') +
    theme(axis.text.x = element_blank()) +
    labs(y='Cummulative percentage (%)', x = '') + 
    theme(axis.title = element_text(size = 60), legend.title = element_text(size = 50),
          legend.text = element_text(size=30)) + 
    geom_text(mapping = aes(x= gap_pattern, y=csum+250, label = csum_percent), angle=90 )
  outpdf(infile, 'motif_csum_intra', width = 32, height = 15)
  
  
  
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
  divdf = data.frame('divs' = c(chao_abinext$Estimator, chao_aginext$Estimator, chao_abinext$Observed, chao_aginext$Observed), 'type' = c('Potential', 'Potential', 'Observed', 'Observed'), 'source' = c('Paratope', 'Epitope', 'Paratope', 'Epitope'))
  print(divdf)
  divdf$source = factor(c('PPI motif', 'PPI motif partner'), levels = c('PPI motif', 'PPI motif partner'))
  ytext = divdf$divs+(0.08*divdf$divs)
  ytext = round(ytext)
  ggplot(data = divdf) +
    geom_bar(stat = 'identity', mapping = aes(reorder(type,divs), divs, fill = source), position = 'dodge') + 
    scale_fill_manual(values = c(ppi_mcolor, ppi_pcolor), name = "Motif source", labels = c('PPI motif', 'PPI motif partner')) + 
    facet_wrap(~source, scales='free') + 
    geom_text(mapping = aes(x=type, y=ytext, label=round(divs))) +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          strip.text = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 20)) + 
    labs(x='', y= 'Interaction motif diversity')
  outpdf(infile, 'chao')
  
}


#run stuff
# ppi_abdb_venn()
# ppi_motif_len()
# ppi_max_gap()
ppi_residue_distribution()
# ppi_residue_correlation()
# ppi_motif_distribution_interdomain()
# ppi_motif_distribution_intradomain()
# chao_diversity()