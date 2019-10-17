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
library('psych')
library('circlize')

theme_set(themeakbar())
my_spectral <- colorRampPalette(brewer.pal(8,'Spectral'))(14)
my_spectral20 <- colorRampPalette(brewer.pal(8,'Spectral'))(20)
print(my_spectral20)
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

sequential_dependency_net = function(infile) {
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  df = read_csv(infile)
  print(dim(df))
  print(df)
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
random_motif_corr('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
random_motif_corr('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# seq_dep_circos('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')
# seq_dep_circos('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv')