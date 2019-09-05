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

ppi_mcolor_rgb = c(col2rgb('orange'))/255
ppi_mcolor = rgb(ppi_mcolor_rgb[1], ppi_mcolor_rgb[2], ppi_mcolor_rgb[3])
ppi_pcolor_rgb = c(col2rgb('darkorange'))/255
ppi_pcolor = rgb(ppi_pcolor_rgb[1], ppi_pcolor_rgb[2], ppi_pcolor_rgb[3])



sequential_dependency_net = function(infile) {
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_XXX_X1X_X2X_XX_edge_next.csv'
  df = read_csv(infile)
  print(dim(df))
  print(df)
  node_size_scaler = 1 
  if (grepl('motifpartner', infile)){
    node_color = ppi_pcolor
    node_size_scaler = 1
  }
  else{
    node_color = ppi_mcolor
  }
  motifs = c('XXX', 'XX', 'X1X', 'X2X')
  for (item in motifs) {
    motif = item
    print(motif)
    dfm = df[df$motif == motif,]
    print(dim(dfm))
    name_parts = paste(strsplit(infile, '_')[[1]][5:10], collapse = '_')
    mtag = strsplit(infile, '_')[[1]][15]
    print(mtag)
    outname = sprintf('abdb_figures_2019/%s_%s_%s_seqnet.pdf', name_parts,mtag, motif)
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
    nodes = data.frame(id=st_residue_count$residues, score = as.integer(st_residue_count$n_norm), stringsAsFactors=FALSE, group=st_residue_count$motif)
    edges = data.frame(source= source_target_count$source, target= source_target_count$target, weight = source_target_count$n_norm, stringsAsFactors=FALSE)
    # print(head(nodes))
    print(head(edges))
    edges_out = edges %>% mutate(motif_source = rep(paste0(mtag, '_', motif), nrow(edges))) %>% mutate(source_target = paste0(source, '_', target))
    print(head(edges_out))
    print(outname)
    outnamedf = sprintf('abdb_outfiles_2019/%s_%s_%s_seqnet_edges.csv', name_parts,mtag, motif)
    write_csv(edges_out, outnamedf)
    # stop()
    deleteAllNetworks()
    print(getNodeShapes())
    setNodeShapeDefault('ELLIPSE')
    createNetworkFromDataFrames(nodes, edges)
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

## norm function
normalize_my <- function(x)
{
  return((x- min(x)) /(max(x)-min(x)))
}

seqnet_clustering = function(){
  infiles = list.files('abdb_outfiles_2019', pattern = 'seqnet_edges', full.names = TRUE)
  infiles = infiles[!grepl('homo|mus',infiles)]
  df = infiles %>%map_df(~fread(.)) 
  print(df)
  dfscaled = df %>% group_by(motif_source) %>% mutate(weight = normalize_my(weight))
  print(df[df$weight>1])
  print(dfscaled)
  sdf = spread(dfscaled[,3:5], key = source_target, value = weight)
  # create row annotatation
  print(rownames(sdf))
  metarow = data.frame(
    Motif = c(rep('Epitope',4), rep('Paratope',4), rep('PPI motif',4), rep('PPI motif partner',4)), row.names = sdf$motif_source
  )
  metrowcol = list(Motif = c(Epitope=agcolor, Paratope=abcolor3, 'PPI motif'=ppi_mcolor, 'PPI motif partner' = ppi_pcolor))
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
  outname = 'abdb_figures_2019/paratope_epitope_ppi_seqnet_edges_heatmap.pdf'
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
                  legend_labels =  c('0.2', '0.4', '0.6', '0.8','Sequential\nfrequency\n1'),
                  clustering_distance_rows = 'correlation'
                  )
  pdf(outname, width = 200, height = 30)
  grid.draw(hmap$gtable)
  dev.off()
  system(sprintf('open %s', outname))
  
}


edges_overlap = function(infile) {
  df = read_csv(infile)
  print(df)
  df$origin = sprintf('epitope', seq(1, dim(df)[1])) 
  df$source_target = paste0(df$source, '_',df$target)
  df$motif_origin = paste0(df$motif, '_',df$origin)
  print(df)
  if (grepl('ppimotif', infile)){
    filcol = ppi_mcolor
    mtag = 'ppimotif'
    print(mtag)
  }
  else if (grepl('ppimotifpartner', infile)){
  filcol = ppi_pcolor
  mtag  = 'ppimotifpartner'
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


edges_overlap_motifwise = function() {
  infiles = list.files('abdb_outfiles_2019', pattern = 'edge_next.csv', full.names = TRUE)
  infiles = infiles[grepl('ppi', infiles)]
  print(infiles)
  mdf = tibble()
  hmcolors = c(ppimotif = ppi_mcolor, ppimotifpartner = ppi_pcolor)
  for (infile in infiles){
    parts = strsplit(infile, split = '_')[[1]]
    print(parts)
    source_type = parts[15]
    print(source_type)
    df = read_csv(infile)
    df$source_type = rep(source_type, nrow(df))
    mdf = rbind(mdf, df)
  }
  print(mdf)
  motifs = unique(mdf$motif)
  for (motif in motifs){
    print(motifs)
    df = mdf[mdf$motif == motif,]
    df$source_target = paste0(df$source, '_',df$target)
    df$motif_origin = paste0(df$motif, '_',df$source_type)
    motifs = unique(df$motif_origin)
    sources = unique(df$source_type)
    print(sources)
    x = sapply(motifs, function(x) subset(df, motif_origin == x)$source_target)
    label_count = sprintf('%s: %s', sources, sapply(sources, function(x) length(unique(subset(df, source_type==x)$source_target))))
    print(length(x))
    print(label_count)
    nproteindf = proteins_in_motif()
    print(nproteindf)
    nprotein1 = nproteindf[nproteindf$gap_pattern1==motif,]$nproteins
    nprotein2 = nproteindf[nproteindf$gap_pattern1==motif,]$nproteins2
    sum_label_count = sprintf('Total edges: %s; Proteins: %s, %s', length(unique(df$source_target)), nprotein1, nprotein2)
    print(sum_label_count)
    # stop()
    parts = strsplit(infile, '_')[[1]]
    nparts = length(parts)
    print(nparts)
    outname = sprintf('abdb_figures_2019/ppimotif_ppimotifpartner_motifwise_seqdep_venn_%s.tiff', motif)
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
                              # alpha=c(0.9,0.7, 0.5, 0.2),
                              alpha=c(0.9,0.7),
                              # rotation.degree = 180,
                              cex = 5, 
                              cat.fontface=1,
                              cat.cex = 3,
                              cat.dist = c(-0.1,-0.1),
                              # cat.pos = c(0.5, 1),
                              # cat.just = list(c(0.2, 10), c(0.2,10)),
                              cat.col = agcolor2,
                              category.names = label_count,
                              main = sum_label_count, 
                              main.col = agcolor2,
                              main.pos = c(0.5,0.05),
                              main.fontface = 1,
                              main.cex = 3
                              )
    system(sprintf('open %s', outname))
  }
}

proteins_in_motif = function(){
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv'
  inmotif = c('XXX', 'XX', 'X1X', 'X2X')
  df0 = read_csv(infile)
  print(df0)
  df = read.csv(infile) %>% group_by(gap_pattern1) %>% summarise(nproteins = length(pdbid1))
  df  = df[df$gap_pattern1 %in% inmotif,]
  print(df)
  df2 = read.csv(infile) %>% group_by(gap_pattern2) %>% summarise(nproteins = length(pdbid2))
  df2  = df2[df2$gap_pattern2 %in% inmotif,]
  print(df2)
  df$nproteins2 = df2$nproteins
  return = df
}

# run stuff
# sequential_dependency_net('abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_XXX_X1X_X2X_XX_ppimotif_edge_next.csv')
# sequential_dependency_net('abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_XXX_X1X_X2X_XX_ppimotifpartner_edge_next.csv')
# seqnet_clustering()
# edges_overlap('abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_XXX_X1X_X2X_XX_ppimotif_edge_next.csv')
# edges_overlap('abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_XXX_X1X_X2X_XX_ppimotifpartner_edge_next.csv')
edges_overlap_motifwise()
# proteins_in_motif()
