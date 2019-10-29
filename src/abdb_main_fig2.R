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


motif.length.distribution_epitope <- function(){
  infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_len.csv'
  df = read.csv(infile)
  print(head(df))
  print(group_by(df, segment))
  df$motif_len = as.numeric(df$motif_len)
  meds = df %>% group_by(segment) %>% summarise(med = median(motif_len))
  print(meds)
  med_labels = paste0('Median: ', meds$med)
  print(med_labels)
  # stop()
  ggplot(data=df) + 
    # geom_bar(mapping = aes(x=motif_len), fill=my_spectral[1]) + 
    # geom_bar(stat='count', mapping = aes(x=motif_len), fill='red') + 
    geom_density(mapping = aes(x=motif_len), fill=alpha(agcolor,0.99), color=alpha(agcolor, 0.99)) + 
    facet_wrap(~ segment, ncol = 3, scales = 'free') + 
    labs(x='Interaction motif length', y='Density') + 
    geom_text(data = meds, mapping = aes(x=med, y=0, label=med_labels), hjust=0) +
    xlim(0,10) +
    coord_flip()
  outname = 'abdb_figures_2019/motif_length_distribution_epitope.pdf'
  ggsave(outname, height = 9, width=7)
  system(sprintf('open %s', outname))
}


motif.length.distribution_paratope <- function(){
  infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_len.csv'
  df = read.csv(infile)
  print(head(df))
  meds = df %>% group_by(segment) %>% summarise(med = median(motif_len))
  print(meds)
  med_labels = paste0('Median: ', meds$med)
  print(med_labels)
  ggplot(data=df) + 
    # geom_bar(mapping = aes(x=motif_len), fill=my_spectral[1]) + 
    # geom_bar(stat='count', mapping = aes(x=motif_len), fill='red') + 
    geom_density(mapping = aes(x=motif_len), fill=alpha(abcolor3,0.99), color=alpha(abcolor3, 0.99)) + 
    facet_wrap(~ segment, ncol = 3, scales = 'free') + 
    labs(x='Interaction motif length', y='Density') + 
    geom_text(data = meds, mapping = aes(x=med, y=0, label=med_labels), hjust=0) +
    xlim(0,10) +
    coord_flip()
  outname = 'abdb_figures_2019/motif_length_distribution.pdf'
  ggsave(outname, height = 9, width=7)
  system(sprintf('open %s', outname))
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
    labs(x='Interaction motif length (paratope)', y='Interaction motif length (epitope)')
  outname = 'abdb_figures/ab_ag_motiflen_correlation.pdf'
  ggsave(outname, height = 9, width=6)
  system(sprintf('open %s', outname))
  
}

ab.ag.venn = function(){
  infile_ab = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx.csv' 
  infile_ag = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx.csv' 
  # infile_ag = 'abdb_outfiles/respairs_epitope_segment_notationx_top-1_gap_patterns.csv' 
  dfab = read_csv(infile_ab)
  dfag = read_csv(infile_ag)
  motif_ab = dfab$gap_patterns
  # print(motif_ab)
  motif_ag = dfag$gap_patterns
  x = list(motif_ab,motif_ag)
  partitions = get.venn.partitions(x)
  overlap = partitions$..values..$'1'
  pdf = dfab[dfab$gap_patterns %in% overlap,]
  countpdf = count(pdf, gapstrstatus)
  countpdf$percent = round(countpdf$n/sum(countpdf$n),2)
  print(countpdf)
  stop()
  #random sampling
  # motif_ab_random = sample(motif_ab, 391)
  # print(motif_ab_random)
  # sapply(1:10, function(xet) get.venn.partitions(x))
  ###
  # outname ='abdb_figures/ab_ag_venn_random.tiff'
  inname = strsplit(tail(strsplit(infile_ab, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_venn.tiff', outfigdir,inname)
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

motif_coverage = function(){
  infile = 'abdb_outfiles_2019/motif_coverage_skip30.csv' 
  df = read.csv(infile, stringsAsFactors = FALSE)
  df = transform(df, unique_motif = as.factor(unique_motif))
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
  outname = 'abdb_figures_2019/motif_coverage.pdf'
  ggsave(outname, height = 10, width = 15)
  system(sprintf('open %s', outname))
}

chao_diversity = function(){
  #compute diversity estimate using chao2
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged.csv'
  df = read_csv(infile)
  print(df)
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
  divdf$source = factor(divdf$source, levels = c('Paratope', 'Epitope'))
  ytext = divdf$divs+(0.08*divdf$divs)
  ytext = round(ytext)
  ggplot(data = divdf) +
    geom_bar(stat = 'identity', mapping = aes(reorder(type,divs), divs, fill = source), position = 'dodge') + 
    scale_fill_manual(values = c(abcolor3, agcolor), name = "Motif source", labels = c('Paratope', 'Epitope')) + 
    facet_wrap(~source, scales='free') + 
    geom_text(mapping = aes(x=type, y=ytext, label=round(divs))) +
    labs(x='', y= 'Interaction motif diversity')
  outname = 'abdb_figures_2019/chao_diversity.pdf'
  ggsave(outname, height = 9, width = 7)
  system(sprintf('open %s', outname))
  
}

abdb_resolution = function() {
  infile = 'abdb_outfiles_2019/abdb_resolution.csv'
  df = read_csv(infile)
  # df = subset(df, df$resolution >5)
  # print(df)
  # stop()
  ggplot(data= df) + 
    geom_density(mapping = aes(x=resolution), fill= abcolor3, color = abcolor3) + 
    scale_x_continuous(breaks = seq(0,30, by=2)) + 
    labs(y= 'Density', x= 'Resolution (Angstrom)')
  outname = 'abdb_figures_2019/abdb_resolution_distribution.pdf'
  ggsave(outname)
  system(sprintf('open %s', outname))
}


### merged from old main figure 4
sequential_dependency_net = function(infile) {
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  df = read_csv(infile)
  df$source = df$source2
  df$target = df$target2
  print(dim(df))
  print(df)
  # stop()
  node_size_scaler = 1 
  mtag = 'paratope'
  if (grepl('epitope', infile)){
    node_color = agcolor
    node_size_scaler = 2
    mtag = 'epitope'
  }
  else{
    node_color = abcolor3
  }
  motifs = c('XXX', 'XX', 'X1X', 'X2X')
  for (item in motifs) {
    motif = item
    print(motif)
    dfm = df[df$motif == motif,]
    print(dim(dfm))
    name_parts = paste(strsplit(infile, '_')[[1]][5:10], collapse = '_')
    outname = sprintf('abdb_figures_2019/%s_%s.pdf', name_parts, motif)
    print(outname)
    if(file.exists(outname)){
      response = file.remove(outname)
    }
    cytoscapePing()
    cytoscapeVersionInfo()
    print(getLayoutNames())
    source_count = count(dfm, source, motif)
    target_count = count(dfm, target, motif)
    source_target_count = count(dfm, source, target) %>% mutate(n_norm = (n-min(n))/(max(n)-min(n))*10) # add min max scalling, scale for visibility
    # print(source_count)
    # print(target_count)
    # print(source_target_count)
    source_target = data_frame(residues = c(dfm$source, dfm$target), motif =c(dfm$motif,dfm$motif))
    st_residue_count = count(source_target, residues,motif) %>% mutate(n_norm = (n-min(n)+1)/(max(n)-min(n))*100 ) # add min max scalling, scale for visibility
    print(st_residue_count)
    st_residue_count[st_residue_count <=1] = 1 # replace 0s with 1s, cytoscape node size must not be 0
    print(st_residue_count)
    print(source_target_count)
    # stop()
    print(head(st_residue_count))
    # nodes = data_frame(id = c(source_count$source, target_count$target), group = c(source_count$motif, target_count$motif),
    #                    score = as.integer(c(source_count$n, target_count$n)), stringsAsFactors = FALSE)
    nodes2 = data.frame(id=st_residue_count$residues, score = as.integer(st_residue_count$n_norm), stringsAsFactors=FALSE, group=st_residue_count$motif)
    edges = data.frame(source= source_target_count$source, target= source_target_count$target, weight = source_target_count$n_norm, stringsAsFactors=FALSE)
    # print(head(nodes))
    print(head(edges))
    edges_out = edges %>% mutate(motif_source = rep(paste0(mtag, '_', motif), nrow(edges))) %>% mutate(source_target = paste0(source, '_', target))
    print(head(edges_out))
    print(outname)
    outnamedf = sprintf('abdb_outfiles_2019/%s_%s_%s_seqnet_edges.csv', name_parts,mtag, motif)
    print(outnamedf)
    write_csv(edges_out, outnamedf)
    # stop()
    deleteAllNetworks()
    print(getNodeShapes())
    setNodeShapeDefault('ELLIPSE')
    createNetworkFromDataFrames(nodes2, edges)
    layoutNetwork('attribute-circle')
    print(getLayoutNames())
    node_attribute = getTableColumns()
    # print(node_attribute)
    setNodeSizeBypass(node.names = node_attribute$name, new.sizes = as.numeric(node_attribute$score))
    setEdgeLineWidthMapping('weight', mapping.type = 'p')
    setEdgeTargetArrowShapeDefault('ARROW')
    setNodeColorDefault(node_color)
    setNodeFontSizeDefault(30)
    exportImage(outname, type = 'PDF')
    system(sprintf('open %s', outname))
    deleteAllNetworks()
  }
}

sequential_dependency_cluster = function(){
  infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  infile2 = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  dfag = read_csv(infile)
  dfag$origin = sprintf('epitope', seq(1, dim(dfag)[1])) 
  dfag$source_target = paste0(dfag$source, '_',dfag$target)
  dfag$motif_origin = paste0(dfag$motif, '_',dfag$origin)
  dfab = read_csv(infile2)
  dfab$origin = sprintf('paratope', seq(1, dim(dfab)[1])) 
  dfab$source_target = paste0(dfab$source, '_',dfab$target)
  dfab$motif_origin = paste0(dfab$motif, '_',dfab$origin)
  print(dfab)
  print(dfag)
  ins = c('XXX', 'XX')
  dfab = dfab[dfab$motif %in% ins,]
  dfag = dfag[dfag$motif %in% ins,]
  # abag_inter = intersect(dfab$source_target, dfag$source_target)
  mdf = rbind(dfab,dfag)
  print(dim(mdf))
  motif_origins = unique(mdf$motif_origin)
  motif_intersect = subset(mdf, mdf$motif_origin==motif_origins[1])$source_target
  print(length(motif_intersect))
  for (item in motif_origins){
    source_target = subset(mdf, mdf$motif_origin==item)$source_target
    print(length(source_target))
    temp = intersect(motif_intersect,source_target)
    motif_intersect = temp 
  }
  print(motif_intersect)
  stop()
  # print(dim(mdf))
  # mdf = mdf[mdf$source_target %in% abag_inter,]
  # print(dim(mdf))
  count_mdf = count(mdf, source_target, motif_origin)
  scount_mdf = spread(count_mdf, key = source_target, value = n)
  print(count_mdf)
  print(scount_mdf)
  scount_mdfnona = drop_na(scount_mdf)
  print(scount_mdfnona)
}


sequential_dependency_cluster2 = function(infile){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  df = read_csv(infile)
  # ins  = c('XXX', 'XX', 'X1X')
  # df  = df[df$motif %in% ins,]
  df$origin = sprintf('epitope', seq(1, dim(df)[1])) 
  df$source_target = paste0(df$source2, '_',df$target2)
  df$motif_origin = paste0(df$motif, '_',df$origin)
  countdf = count(df, source_target, motif) %>% group_by(motif) %>% mutate(n = (n-min(n))/(max(n)-min(n)))
  sdf = spread(countdf, key=source_target, value=n)
  print(sdf)
  edge_intersect = colSums(is.na(sdf[,1:ncol(sdf)])) == 0
  print(edge_intersect)
  # print(edge_intersect)
  sdf[is.na(sdf)] = 0
  # sdf2 = sdf
  sdf2 = sdf[, edge_intersect]
  # stop()
  sdf2mat = sdf2[,-1]
  print(sdf2mat)
  name_parts = paste(strsplit(infile, '_')[[1]][5:10], collapse = '_')
  outname = sprintf('abdb_figures_2019/%s_seqdep_heatmap.pdf', name_parts)
  print(sdf2$motif)
  hmap = pheatmap(sdf2mat,
                  labels_row = sdf2$motif,
                  # color = my_spectral[1:14],
                  color = viridis(7),
                  number_color = 'white',
                  cellwidth = 30, cellheight = 30,
                  border_color = 'white',
                  display_numbers = TRUE,
                  cutree_rows = 4,
                  legend = TRUE,
                  legend_breaks = c(0.2,0.4,0.6,0.8,1),
                  legend_labels =  c('0.2', '0.4', '0.6', '0.8','Sequential frequency\n1')
                  )
  pdf(outname, width =20, height = 4)
  grid.draw(hmap$gtable)
  dev.off()
  system(sprintf('open %s', outname))
}


sequential_dependency_cluster2_alledges = function(infile){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  df = read_csv(infile)
  # ins  = c('XXX', 'XX', 'X1X')
  # df  = df[df$motif %in% ins,]
  df$origin = sprintf('epitope', seq(1, dim(df)[1])) 
  df$source_target = paste0(df$source, '_',df$target)
  df$motif_origin = paste0(df$motif, '_',df$origin)
  # countdf = count(df, source_target, motif)
  countdf = count(df, source_target, motif) %>% group_by(motif) %>% mutate(n = (n-min(n))/(max(n)-min(n)))
  sdf = spread(countdf, key=source_target, value=n)
  edge_intersect = colSums(is.na(sdf[,1:ncol(sdf)])) == 0
  # print(edge_intersect)
  sdf[is.na(sdf)] = 0
  sdf2 = sdf
  # sdf2 = sdf[, edge_intersect]
  sdf2mat = sdf2[,-1]
  print(sdf2mat)
  name_parts = paste(strsplit(infile, '_')[[1]][5:10], collapse = '_')
  outname = sprintf('abdb_figures_2019/%s_seqdep_heatmap_alledges.pdf', name_parts)
  print(sdf2$motif)
  hmap = pheatmap(sdf2mat,
                  labels_row = sdf2$motif,
                  # scale = 'row',
                  color = my_spectral[1:14],
                  cellwidth = 10, cellheight = 10,
                  border_color = 'white',
                  # display_numbers = TRUE,
                  cutree_rows = 4,
                  width = 8, height = 4,
                  )
  pdf(outname, width = 45, height = 4)
  grid.draw(hmap$gtable)
  dev.off()
  system(sprintf('open %s', outname))
}

edges_overlap = function(infile) {
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  df = read_csv(infile)
  df$origin = sprintf('epitope', seq(1, dim(df)[1])) 
  df$source_target = paste0(df$source2, '_',df$target2)
  df$motif_origin = paste0(df$motif, '_',df$origin)
  print(df)
  filcol = abcolor3
  if (grepl('epitope', infile)) {filcol = agcolor}
  motifs = unique(df$motif)
  x = sapply(motifs, function(x) subset(df, motif == x)$source_target)
  label_count = sprintf('%s: %s', motifs, sapply(motifs, function(x) length(unique(subset(df, motif==x)$source_target))))
  print(length(x))
  print(label_count)
  sum_label_count = sprintf('Total edges: %s', length(unique(df$source_target)))
  print(sum_label_count)
  # stop()
  name_parts = paste(strsplit(infile, '_')[[1]][5:10], collapse = '_')
  outname = sprintf('abdb_figures_2019/%s_seqdep_venn.tiff', name_parts)
  print(outname)
  venn.plot <- venn.diagram(x=x, filename=outname,
                            output = FALSE ,
                            imagetype="tiff" ,
                            # height = 800 , 
                            # width = 900 , 
                            resolution = 300,
                            compression = "lzw",
                            lwd = 2,
                            lty = 'blank',
                            fill = c(filcol, filcol, filcol, filcol),
                            # alpha=c(0.2,0.5, 0.7, 0.9),
                            alpha=c(0.9,0.7, 0.5, 0.2),
                            # rotation.degree = 180,
                            cex = 5, 
                            cat.fontface=1,
                            cat.cex = 3,
                            # cat.dist = c(-0.05,-0.15),
                            cat.just = list(c(0, 21), c(1,21), c(0.9,1), c(0,1)),
                            cat.col = c(filcol, filcol, filcol, filcol),
                            category.names = label_count,
                            main = sum_label_count, 
                            main.col = filcol,
                            main.pos = c(0.5,0.05),
                            main.fontface = 1,
                            main.cex = 3
                            )
  system(sprintf('open %s', outname))
  
}

random_motif_corr = function(infile){
  # correlate observed pair frequency and randomly drawn pair frequence
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  df = read_csv(infile)
  motifs = unique(df$motif)
  print(motifs)
  outdf = data_frame()
  ksdf =data_frame()
  for (motif in motifs){
    mdf = df[df$motif == motif,] #%>% group_by(source, target) %>% count()
    mdf$pairs = paste0(mdf$source2, '-', mdf$target2)
    rdf = data_frame()
    n = dim(mdf)[1]
    nrep = 100
    for (i in seq(1,nrep)){
      nrdf = sample_n(df, n)
      nrdf$rep = rep(i, dim(nrdf)[1])
      rdf = rbind(rdf, nrdf)
    }
    rdf$pairs = paste0(rdf$source, '-', rdf$target)
    countrdf = count(rdf, pairs, rep)
    # countrdf$rep2= countrdf$rep
    # grdf =  spread(countrdf, key = rep, value = n)
    # print(grdf)
    # ct = psych::corr.test(grdf[,2:11])$ci
    # print(ct)
    rdf = rdf[rdf$pairs %in% mdf$pairs,]
    countmdf = count(mdf, pairs)
    countmdf$rep = rep(0, dim(countmdf)[1])
    countrdf = count(rdf, pairs, rep)
    mergeddf = rbind(countrdf, countmdf) %>% spread(key = rep, value = n)
    print(dim(mergeddf))
    # stop()
    ct = psych::corr.test(mergeddf[2:dim(mergeddf)[2]])$ci
    ct$motif = rep(motif, dim(ct)[1])
    ct$test_type = c(rep('ovr',nrep), rep('rvr', dim(ct)[1]-nrep))
    print(ct)
    ks = ks.test(ct[ct$test_type == 'ovr',]$r, ct[ct$test_type == 'rvr',]$r)
    ksvect = data_frame(motif = motif, dval = ks$statistic, pval = ks$p.value)
    ksdf = rbind(ksvect, ksdf) 
    print(ksvect)
    # stop()
    outdf = rbind(outdf, ct)
    print(ct)
    # stop()
    print(dim(ct))
  }
  print(head(outdf))
  filcols= c(abcolor3,  abcolor3)
  if (grepl('epitope', infile)){
    filcols = c('black', agcolor) 
  }
  ksdf[ksdf == 0] =  '<10e-11'
  print(ksdf)
  # stop()
  ggplot(data=outdf) + 
    # geom_point(mapping = aes(x = n.x, y= n.y)) +
    geom_density(data = outdf, mapping = aes(r, fill=test_type), color =NA) + 
    facet_wrap(~motif, scales = 'free') + 
    labs(x= 'Pearson cor.', y = 'Density') + 
    scale_fill_manual(values = c(filcols[1], ggplot2::alpha(filcols[2],0.4)),
                      labels = c('Observed vs random', 'Random vs random')) +
    theme(legend.position = 'bottom',
          axis.title = element_text(size=30),
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          strip.text = element_text(size = 20),
          axis.text = element_text(size = 20)) + 
    geom_text(data=ksdf, mapping = aes(x=-0.5, y=-1, label = sprintf('KS test p-val: %s',formatC(pval,digits = 3))), size = 6, hjust=0, vjust=0)
  outpdf(infile, 'density')
}

seq_dep_circos = function(infile){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  df = read_csv(infile)
  df2 = read_csv('abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil_pc.csv')
  print(df2)
  motifs = c('XXX', 'XX', 'X1X', 'X2X')
  seqdf = df2 %>% group_by(ab_motif) %>% summarise(len = length(paratope)) 
  print(seqdf[order(seqdf$len, decreasing = TRUE ),])
  seqdf2 = df2 %>% group_by(ag_motif) %>% summarise(len = length(paratope)) 
  print(seqdf2[order(seqdf2$len, decreasing = TRUE ),])
  stop()
  print(df)
  countdf = count(df, source2, target2, motif)
  print(countdf)
  aa = unique(countdf$source2)
  aa2 = c(aa, paste0(aa,"'"))
  print(aa)
  grid.col = c(my_spectral20, my_spectral20)
  names(grid.col) = aa2
  motifs = unique(countdf$motif)
  for (motif in motifs){
    motifdf = countdf[countdf$motif== motif,]
    # motifdf = motifdf[motifdf$n>5,]
    print(motifdf)
    # stop()
    inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
    outname = sprintf('%s/%s_%s_%s.pdf',outfigdir, inname, motif,'circos')
    pdf(outname)
    chordDiagram(motifdf[,c(1,2,4)], grid.col = grid.col,
                 scale=FALSE,
                 annotationTrack = c('grid'),
                 annotationTrackHeight = c(0.15),)
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1]+0.5, CELL_META$sector.index, 
                  cex = 2)
    }, bg.border = NA)  
    system(sprintf('open %s', outname))
    circos.clear()
    dev.off()
    # stop()
  }
}


# runstuff
# motif.length.distribution_epitope()
# motif.length.distribution_paratope()
# ab_ag_motiflen_correlation()
ab.ag.venn()
# motif_coverage()
# chao_diversity()
# abdb_resolution()


# run stuff
# note to self: cytoscape is finicky: need to run 2x, separately, in order for the program to update the plots to new changes.
# sequential_dependency_net('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# sequential_dependency_net('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# sequential_dependency_cluster()
# sequential_dependency_cluster2('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# sequential_dependency_cluster2_alledges('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# sequential_dependency_cluster2('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# sequential_dependency_cluster2_alledges('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# edges_overlap('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# edges_overlap('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# random_motif_corr('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# random_motif_corr('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# seq_dep_circos('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# seq_dep_circos('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
  