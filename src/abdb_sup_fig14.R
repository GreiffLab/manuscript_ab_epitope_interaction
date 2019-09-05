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
greys = brewer.pal(10,'Greys')
print(greys)
agcolor2 = greys[2]
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
  print(df)
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


seqnet_clustering = function(){
  infiles = list.files('abdb_outfiles_2019', pattern = 'seqnet_edges', full.names = TRUE)
  infiles = infiles[grepl('homo|mus',infiles)]
  df = infiles %>%map_df(~fread(.)) 
  print(unique(df$motif_source))
  dfscaled = df %>% group_by(motif_source) %>% mutate(weight = normalize_my(weight))
  print(df[df$weight>1])
  print(dfscaled)
  sdf = spread(dfscaled[,3:5], key = source_target, value = weight)
  print(sdf)
  stop()
  # create row annotatation
  print(rownames(sdf))
  metarow = data.frame(
    Motif = c(rep('epitope',4), rep('paratope',4), rep('ppimotif',4), rep('ppimotifpartner',4)), row.names = sdf$motif_source
  )
  metrowcol = list(Motif = c(epitope=agcolor, paratope=abcolor3, ppimotif=ppi_mcolor, ppimotifpartner = ppi_pcolor))
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
  hmap = pheatmap(sdf2mat,
                  labels_row = sdf$motif_source,
                  annotation_row = metarow,
                  annotation_colors = metrowcol[1],
                  color = viridis(7),
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
  grid.draw(hmap$gtable)
  dev.off()
  system(sprintf('open %s', outname))
  
}



# run stuff
# sequential_dependency_net_homo_mus('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_homo_XXX_X1X_X2X_XX_edge_next.csv')
# sequential_dependency_net_homo_mus('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_mus_XXX_X1X_X2X_XX_edge_next.csv')
# sequential_dependency_net_homo_mus('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_homo_XXX_X1X_X2X_XX_edge_next.csv')
# sequential_dependency_net_homo_mus('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_mus_XXX_X1X_X2X_XX_edge_next.csv')
seqnet_clustering()




