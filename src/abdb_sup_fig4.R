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

human_mouse_overlap = function(){
  infile  = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species.csv'
  df = read_csv(infile)
  print(df)
  print(unique(df$amolecule))
  mousedf =  subset(df,hspecies == "MUS MUSCULUS" )
  humandf =  subset(df,hspecies == "HOMO SAPIENS" )
  x= list('Paratope mouse' = mousedf$ab_motif, "Epitope mouse" = mousedf$ag_motif, 'Paratope human' = humandf$ab_motif, 'Epitope human' = humandf$ag_motif)
  print(x$ab_human)
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_venn.tiff', outfigdir,inname)
  print(outname)
  print(names(x))
  # stop()
  venn.plot <- venn.diagram(x=x, filename=outname,
                            output = FALSE ,
                            imagetype="tiff" ,
                            # height = 800 , 
                            # width = 900 , 
                            resolution = 300,
                            compression = "lzw",
                            lwd = 2,
                            lty = 'blank',
                            fill = c(abcolor3, agcolor, abcolor3, agcolor),
                            alpha=c(0.7,0.2, 0.7, 0.3),
                            # rotation.degree = 180,
                            cex = 5, 
                            cat.fontface=3,
                            cat.cex = 3,
                            # cat.dist = c(-0.05,-0.15),
                            cat.just = list(c(0, 21), c(1,21), c(0.9,1), c(0,1)),
                            cat.col = c(abcolor3, 'grey', abcolor3, 'grey'),
                            # category.names = names(x)
                            )
  system(sprintf('open %s', outname))
  
}


paratope_top_antigen = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
  df = read_csv(infile)
  print(head(df))
  n=20
  top_nmotif = count(df, ab_motif, sort = TRUE) %>% top_n(n)
  print(top_nmotif)
  df = df[df$ab_motif %in% top_nmotif$ab_motif,]
  print(df)
  countdf = count(df, ab_motif, ag_species2)
  print(head(countdf))
  countdf2 = count(countdf, ab_motif) 
  toptendf = countdf2[order(-countdf2$nn),][1:20,]
  print(toptendf)
  labels = paste0('',toptendf$nn)
  print(labels)
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_shared_paratope_motif_antigen.pdf', outfigdir,inname)
    ggplot(data = toptendf) + 
    geom_bar(mapping = aes(x=reorder(ab_motif,-nn), y = nn), stat = 'identity', fill=abcolor3) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    geom_text(mapping= aes(x = toptendf$ab_motif, y=toptendf$nn+1, label=labels)) +
    labs(x='Top 20 paratope motif', y='# of antigen classes')
  ggsave(outname, width = 10, height = 8)
  system(sprintf('open %s', outname))
}


paratope_top_complex = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
  df = read_csv(infile)
  print(head(df))
  n=20
  # stop()
  top_nmotif = count(df, ab_motif, sort = TRUE) %>% top_n(n)
  print(top_nmotif)
  df = df[df$ab_motif %in% top_nmotif$ab_motif,]
  print(df)
  countdf = count(df, ab_motif, pdbid, hspecies)
  print(head(countdf))
  countdf2 = count(countdf, ab_motif) 
  toptendf = countdf2[order(-countdf2$n),][1:20,]
  print(toptendf)
  labels = paste0('',toptendf$n)
  print(labels)
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_shared_paratope_motif_complex.pdf', outfigdir,inname)
    ggplot(data = toptendf) + 
    geom_bar(mapping = aes(x=reorder(ab_motif,-n), y = n), stat = 'identity', fill=abcolor3) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    geom_text(mapping= aes(x = toptendf$ab_motif, y=toptendf$n+10, label=labels)) +
    labs(x='Top 20 paratope motif', y='# of structures')
  ggsave(outname, width = 10, height = 8)
  system(sprintf('open %s', outname))
}


epitope_top_complex = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
  df = read_csv(infile)
  print(head(df))
  n=20
  # stop()
  top_nmotif = count(df, ag_motif, sort = TRUE) %>% top_n(n)
  print(top_nmotif)
  df = df[df$ag_motif %in% top_nmotif$ag_motif,]
  print(df)
  countdf = count(df, ag_motif, pdbid)
  print(head(countdf))
  countdf2 = count(countdf, ag_motif) 
  toptendf = countdf2[order(-countdf2$nn),][1:20,]
  print(toptendf)
  labels = paste0('',toptendf$nn)
  print(labels)
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_shared_epitope_motif_complex.pdf', outfigdir,inname)
    ggplot(data = toptendf) + 
    geom_bar(mapping = aes(x=reorder(ag_motif,-nn), y = nn), stat = 'identity', fill=agcolor) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    geom_text(mapping= aes(x = toptendf$ag_motif, y=toptendf$nn+10, label=labels)) +
    labs(x='Top 20 paratope motif', y='# of structures')
  ggsave(outname, width = 10, height = 8)
  system(sprintf('open %s', outname))
}


cdr_fr_overlap_ab = function(infile) {
  # infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2_ab_motif_overlap_min_max.csv'
  df = read_csv(infile)
  mat = as.matrix(df[seq(1,dim(df)[1]),seq(2,dim(df)[2])])
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_heatmap.pdf', outfigdir,inname)
  pheatmap(mat = mat, filename = outname, 
           color = colorRampPalette(brewer.pal(n = 3, name ="Blues"))(100),
           border_color = NA,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           labels_row = df$rownames,
           display_numbers = TRUE,
           number_color = 'Black',
           number_format = '%i')
  system(sprintf('open %s', outname))
}


cdr_fr_overlap_ag = function(infile) {
  # infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2_ab_motif_overlap_min_max.csv'
  df = read_csv(infile)
  mat = as.matrix(df[seq(1,dim(df)[1]),seq(2,dim(df)[2])])
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_heatmap.pdf', outfigdir,inname)
  pheatmap(mat = mat, filename = outname, 
           color = colorRampPalette(brewer.pal(n = 3, name ="Greys"))(100),
           border_color = NA,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           labels_row = df$rownames,
           display_numbers = TRUE,
           number_color = 'Black',
           number_format = '%i')
  system(sprintf('open %s', outname))
}


paratope_top_imgt_vgene = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2_imgt_vgene.csv'
  df = read_csv(infile)
  df[df == "Mus musculus_C57BL/6"] = 'Mus musculus (C57BL/6)'
  n=20
  top_nmotif = count(df, ab_motif, sort = TRUE) %>% top_n(n)
  print(top_nmotif)
  # stop()
  df = df[df$ab_motif %in% top_nmotif$ab_motif,]
  vgene_spec = count(df, imgt_vgenename, imgt_species, abchain)
  vgene_spec_count = count(vgene_spec, imgt_species, abchain)
  motifcount = count(df, ab_motif)
  topmotif = motifcount[order(-motifcount$n),]
  print(topmotif)
  # countdf = count(df, ab_motif, imgt_vgenename, imgt_species,abchain)
  countdf = count(df, ab_motif, imgt_vgenename, imgt_species)
  print(head(countdf))
  countdf2 = count(countdf, ab_motif, imgt_species) 
  # countdf2 = count(countdf, ab_motif, imgt_species,abchain) 
  toptendf = countdf2[order(-countdf2$nn),]
  print(toptendf)
  labels = paste0('',toptendf$nn)
  print(labels)
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_shared_paratope_motif_imgt_vgene.pdf', outfigdir,inname)
    ggplot(data = toptendf) + 
    geom_bar(mapping = aes(x=reorder(ab_motif,-nn), y = nn, fill = imgt_species), stat = 'identity', position = 'dodge') + 
    theme(axis.text.x = element_text(angle = 90)) +
    # facet_wrap(~ abchain) +
    geom_text(mapping= aes(x = toptendf$ab_motif, y=toptendf$nn+1, label=labels, fill=imgt_species), position = position_dodge(width=1)) +
    scale_fill_manual(values = c(alpha(abcolor3, 0.5), abcolor3)) +
    labs(x='Top 20 paratope motif', y='# of unique V genes', fill='Species (IMGT naming)') + 
    # put legend inside
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    ) 
  ggsave(outname, width = 10, height = 8)
  system(sprintf('open %s', outname))
}


paratope_top20_cumulative = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
  df = read_csv(infile)
  print(head(df))
  n=20
  countdf = count(df, ab_motif)
  print(countdf)
  total_motif = sum(countdf$n)
  total_motif_label = sprintf('Total number of motifs: %s', total_motif)
  top_nmotif = count(df, ab_motif, sort = TRUE) %>% top_n(n) %>% mutate(csum = round(cumsum(n)/dim(df)[1]*100, 1))
  print(top_nmotif)
  print(cumsum(top_nmotif$n))
  # stop()
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_paratope_top20_cumulative_percent.pdf', outfigdir,inname)
    ggplot(data = top_nmotif) + 
    geom_bar(mapping = aes(x=reorder(ab_motif,csum), y = csum), stat = 'identity', fill=abcolor3) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    geom_text(mapping= aes(x = ab_motif, y=csum+3, label=csum), angle=90) +
    geom_text(data = df, mapping= aes(x = 1, y= 80, label=total_motif_label), hjust=0) +
    labs(x='Top 20 paratope motif', y='Cumulative percentage (%)')
  ggsave(outname, width = 5, height = 8)
  system(sprintf('open %s', outname))
}


epitope_top20_cumulative = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
  df = read_csv(infile)
  print(head(df))
  n=20
  countdf = count(df, ag_motif)
  print(countdf)
  total_motif = sum(countdf$n)
  total_motif_label = sprintf('Total number of motifs: %s', total_motif)
  top_nmotif = count(df, ag_motif, sort = TRUE) %>% top_n(n) %>% mutate(csum = round(cumsum(n)/dim(df)[1]*100, 1))
  print(top_nmotif)
  print(cumsum(top_nmotif$n))
  # stop()
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_epitope_top20_cumulative_percent.pdf', outfigdir,inname)
    ggplot(data = top_nmotif) + 
    geom_bar(mapping = aes(x=reorder(ag_motif,csum), y = csum), stat = 'identity', fill=agcolor) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    geom_text(mapping= aes(x = ag_motif, y=csum+3, label=csum), angle=90) +
    geom_text(data = df, mapping= aes(x = 1, y= 70, label=total_motif_label), hjust=0) +
    labs(x='Top 20 epitope motif', y='Cumulative percentage (%)')
  ggsave(outname, width = 5, height = 8)
  system(sprintf('open %s', outname))
}


# run stuff
# human_mouse_overlap()
# paratope_top_antigen()
paratope_top_complex()
# cdr_fr_overlap_ab('abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2_ab_motif_overlap_min_max.csv')
# cdr_fr_overlap_ag('abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2_ag_motif_overlap_min_max.csv')
# paratope_top_imgt_vgene()
# epitope_top_complex()
# paratope_top20_cumulative()
# epitope_top20_cumulative()
#overlap equations
#Motif_i \cap Motif_j/\textrm{max}(Motif_i, Motif_j)
#Motif_i \cap Motif_j/\textrm{min}(Motif_i, Motif_j)

