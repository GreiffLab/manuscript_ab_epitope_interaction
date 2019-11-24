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
greys = brewer.pal(9,'Greys')
print(greys)
agcolor2 = greys[7]
# agcolor = greys[4]
blues2 = brewer.pal(6,'Blues')
blues3 = c('#9010FF', '#1010FF', '#C8C8FF')
blues4 = c(blues3, blues2[2:6])
cdrfr_colors = c(oranges,blues)
rewcolors = c(greys[3:5], c('#9010FF', '#1010FF', '#C8C8FF'))
rewcolors2  = c(greys[3:4], blues3[1:2])
blues5 = c('#1010FF','#1010FF','#1010FF',greys[3],'#1010FF')
outfigdir = 'abdb_figures_2019'

ppi_mcolor_rgb = c(col2rgb('orange'))/255
ppi_mcolor = rgb(ppi_mcolor_rgb[1], ppi_mcolor_rgb[2], ppi_mcolor_rgb[3])
ppi_pcolor_rgb = c(col2rgb('darkorange'))/255
ppi_pcolor = rgb(ppi_pcolor_rgb[1], ppi_pcolor_rgb[2], ppi_pcolor_rgb[3])
mouseabcolor = my_spectral[12]

sequential_dependency_net_homo_mus = function(infile) {
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  df = read_csv(infile)
  print(dim(df))
  df$source = df$source2
  df$target = df$target2
  print(df)
  # stop()
  node_size_scaler = 1 
  if (grepl('paratope', infile) & grepl('homo', infile)){
    node_color = abcolor3
    mtag = 'hparatope'
    print(mtag)
  }
  else if (grepl('paratope', infile) & grepl('mus', infile)){
  node_color = mouseabcolor
  mtag  = 'mparatope'
  print(mtag)
  }
  else if (grepl('epitope', infile) & grepl('mus', infile)){
  node_color = agcolor2
  mtag  = 'mepitope'
  print(mtag)
  }
  else if (grepl('epitope', infile) & grepl('homo', infile)){
  node_color = agcolor
  mtag  = 'hepitope'
  print(mtag)
  }
  motifs = c('XXX', 'XX', 'X1X', 'X2X')
  for (item in motifs) {
    motif = item
    print(motif)
    dfm = df[df$motif == motif,]
    print(dim(dfm))
    parts = strsplit(infile, '_')[[1]]
    nparts = length(parts)
    print(parts)
    name_parts = paste(strsplit(infile, '_')[[1]][5:nparts-1], collapse = '_')
    name_parts = paste0(name_parts, '_next')
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


min_max_scale = function(x){
  normx = (x-min(x))/(max(x)-min(x))
  return(normx)
}


seqnet_clustering = function(){
  infiles = list.files('abdb_outfiles_2019', pattern = 'seqnet_edges', full.names = TRUE)
  infiles = infiles[grepl('homo|mus',infiles)]
  print(infiles)
  df = infiles %>%map_df(~fread(.))
  print(df)
  # stop()
  print(unique(df$motif_source))
  dfscaled = df %>% group_by(motif_source) %>% mutate(weight = min_max_scale(weight))
  print(df[df$weight>1])
  print(dfscaled)
  sdf = spread(dfscaled[,3:5], key = source_target, value = weight)
  print(sdf)
  # stop()
  # create row annotatation
  print(rownames(sdf))
  metarow = data.frame(
    Motif = c(rep('Epitope homo',4), rep('Paratope homo',4), rep('Epitope mus',4), rep('Paratope mus',4)), row.names = sdf$motif_source
  )
  metrowcol = list(Motif = c('Epitope homo'=agcolor, 'Paratope homo'=abcolor3, 'Epitope mus'=agcolor2, 'Paratope mus' = mouseabcolor))
  print(metrowcol)
  print(metarow)
  # stop()
  print(dim(df))
  sdf[is.na(sdf)] = 0
  sdf2mat = data.frame(sdf[,-1])
  print(sdf2mat)
  rownames(sdf2mat) = sdf$motif_source
  edge_intersect = colSums(is.na(sdf[,1:ncol(sdf)])) == 0
  # sdf2 = sdf[, edge_intersect]
  outname = 'abdb_figures_2019/homo_mus_paratope_epitope_seqnet_edges_heatmap.pdf'
  outnamepng = 'abdb_figures_2019/homo_mus_paratope_epitope_seqnet_edges_heatmap.png'
  hmap = pheatmap(sdf2mat,
                  labels_row = sdf$motif_source,
                  annotation_row = metarow,
                  annotation_colors = metrowcol[1],
                  # color = viridis(7),
                  color = wes_palette('GrandBudapest2')[2:4],
                  number_color = 'white',
                  cellwidth = 30, cellheight = 30,
                  border_color = 'white',
                  display_numbers = TRUE,
                  cutree_rows = 4,
                  legend = TRUE,
                  legend_breaks = c(0.2,0.4,0.6,0.8,1),
                  legend_labels =  c('0.2', '0.4', '0.6', '0.8','Sequential frequency\n1'),
                  clustering_distance_rows = 'correlation'
                  )
  pdf(outname, width = 200, height = 30)
  png(outnamepng, width = 20000, height = 2000, res = 150)
  grid.draw(hmap$gtable)
  dev.off()
  system(sprintf('open %s', outnamepng))
  
}


seq_dep_circos = function(infile){
  df = read_csv(infile)
  df2 = read_csv('abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil_pc.csv')
  print(df)
  # stop()
  motifs = c('XXX', 'XX', 'X1X', 'X2X')
  seqdf = df2 %>% group_by(ab_motif) %>% summarise(len = length(paratope)) 
  print(seqdf[order(seqdf$len, decreasing = TRUE ),])
  seqdf2 = df2 %>% group_by(ag_motif) %>% summarise(len = length(paratope)) 
  print(seqdf2[order(seqdf2$len, decreasing = TRUE ),])
  print(df)
  countdf = count(df, source2, target2, motif)
  print(countdf)
  # stop()
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
    inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
    outname = sprintf('%s/%s_%s_%s.pdf',outfigdir, inname, motif,'circos')
    print(outname)
    # stop()
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


edges_overlap = function(infile) {
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  df = read_csv(infile)
  df$source = df$source2
  df$target = df$target2
  df$origin = sprintf('epitope', seq(1, dim(df)[1])) 
  df$source_target = paste0(df$source, '_',df$target)
  df$motif_origin = paste0(df$motif, '_',df$origin)
  print(df)
  # stop()
  if (grepl('paratope', infile) & grepl('homo', infile)){
    filcol = abcolor3
    mtag = 'hparatope'
    print(mtag)
  }
  else if (grepl('paratope', infile) & grepl('mus', infile)){
  filcol = mouseabcolor
  mtag  = 'mparatope'
  print(mtag)
  }
  else if (grepl('epitope', infile) & grepl('mus', infile)){
  filcol = agcolor2
  mtag  = 'mepitope'
  print(mtag)
  }
  else if (grepl('epitope', infile) & grepl('homo', infile)){
  filcol = agcolor
  mtag  = 'hepitope'
  print(mtag)
  }
  motifs = unique(df$motif)
  x = sapply(motifs, function(x) subset(df, motif == x)$source_target)
  label_count = sprintf('%s: %s', motifs, sapply(motifs, function(x) length(unique(subset(df, motif==x)$source_target))))
  print(length(x))
  print(label_count)
  sum_label_count = sprintf('Total edges: %s', length(unique(df$source_target)))
  print(sum_label_count)
  parts = strsplit(infile, '_')[[1]]
  nparts = length(parts)
  print(nparts)
  name_parts = paste(parts[5:nparts-1], collapse = '_')
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



get_pdbs = function(infile, motif){
 df = read_csv(infile)
 dfm = df[df$gap_patterns == motif,]
 # print(infile)
 parts = strsplit(infile,split = '_')[[1]]
 mtag = paste('PR_',substr(parts[10],1,1))
 stag = strsplit(parts[13], split = '\\.')[[1]][1]
 npdb = paste0(mtag, '_', stag, ': ', nrow(dfm))
 return = npdb
}



edges_overlap_motifwise = function() {
  infiles = list.files('abdb_outfiles_2019', pattern = 'edge_next.csv', full.names = TRUE)
  infiles = infiles[grepl('homo|mus', infiles)]
  infiles2 = list.files('abdb_outfiles_2019', pattern= '_homo.csv|_mus.csv', full.names = TRUE)
  mdf = tibble()
  hmcolors = c(epitope_homo = agcolor, epitope_mus = agcolor2, paratope_homo = abcolor3, paratope_mus = mouseabcolor)
  for (infile in infiles){
    parts = strsplit(infile, split = '_')[[1]]
    print(parts)
    source_type = paste0(parts[10], '_', parts[13])
    print(source_type)
    df = read_csv(infile)
    print(df)
    df$source_type = rep(source_type, nrow(df))
    mdf = rbind(mdf, df)
  }
  print(mdf)
  motifs = unique(mdf$motif)
  for (motif in motifs){
    print(motifs)
    pdbs = sapply(infiles2, function(x) get_pdbs(x,motif))
    print(paste(pdbs[1]))
    # stop()
    df = mdf[mdf$motif == motif,]
    print(df)
    df$source = df$source2
    df$target = df$target2
    # stop()
    df$source_target = paste0(df$source, '_',df$target)
    df$motif_origin = paste0(df$motif, '_',df$source_type)
    motifs = unique(df$motif_origin)
    sources = unique(df$source_type)
    print(sources)
    x = sapply(motifs, function(x) subset(df, motif_origin == x)$source_target)
    label_count = sprintf('%s: %s', sources, sapply(sources, function(x) length(unique(subset(df, source_type==x)$source_target))))
    print(length(x))
    print(label_count)
    sum_label_count = paste0('Total edges: ', 
                              length(unique(df$source_target)), '\n',
                              pdbs[1], '; ', pdbs[2], '; ', 
                              pdbs[3], '; ', pdbs[4])
    print(sum_label_count)
    # stop()
    parts = strsplit(infile, '_')[[1]]
    nparts = length(parts)
    print(nparts)
    outname = sprintf('abdb_figures_2019/homo_mus_motifwise_seqdep_venn_%s.tiff', motif)
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
                              fill = hmcolors,
                              # alpha=c(0.2,0.5, 0.7, 0.9),
                              alpha=c(0.9,0.7, 0.5, 0.2),
                              # rotation.degree = 180,
                              cex = 5, 
                              cat.fontface=1,
                              cat.cex = 3,
                              # cat.dist = c(-0.05,-0.15),
                              cat.just = list(c(0, 21.4), c(1,21.4), c(0.7,1.1), c(0.3,1.1)),
                              cat.col = hmcolors,
                              category.names = label_count,
                              main = sum_label_count, 
                              main.col = agcolor2,
                              main.pos = c(0.5,0.06),
                              main.fontface = 1,
                              main.cex = 2
                              )
    system(sprintf('open %s', outname))
    # stop()
  }
}




# run stuff
# sequential_dependency_net_homo_mus('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_homo_XXX_X1X_X2X_XX_edge_next.csv')
# sequential_dependency_net_homo_mus('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_mus_XXX_X1X_X2X_XX_edge_next.csv')
# sequential_dependency_net_homo_mus('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_homo_XXX_X1X_X2X_XX_edge_next.csv')
# sequential_dependency_net_homo_mus('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_mus_XXX_X1X_X2X_XX_edge_next.csv')
seqnet_clustering()
# seq_dep_circos('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_homo_XXX_X1X_X2X_XX_edge_next.csv')
# seq_dep_circos('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_mus_XXX_X1X_X2X_XX_edge_next.csv')
# seq_dep_circos('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_homo_XXX_X1X_X2X_XX_edge_next.csv')
# seq_dep_circos('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_mus_XXX_X1X_X2X_XX_edge_next.csv')
# edges_overlap('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_homo_XXX_X1X_X2X_XX_edge_next.csv')
# edges_overlap('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_mus_XXX_X1X_X2X_XX_edge_next.csv')
# edges_overlap('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_mus_XXX_X1X_X2X_XX_edge_next.csv')
# edges_overlap('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_homo_XXX_X1X_X2X_XX_edge_next.csv')
# edges_overlap_motifwise()



