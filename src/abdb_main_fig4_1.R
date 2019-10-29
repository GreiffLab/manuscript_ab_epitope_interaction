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
# detach('package:psych')


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
ragcolor_rgb = c(col2rgb('lightsalmon'))/255
ragcolor = rgb(ragcolor_rgb[1], ragcolor_rgb[2], ragcolor_rgb[3])
rabcolor_rgb = c(col2rgb('lightpink'))/255
rabcolor = rgb(rabcolor_rgb[1], rabcolor_rgb[2], rabcolor_rgb[3])


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


paratope_epitope_motif_interaction = function(resize=FALSE){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle.csv'
  name_parts = paste(strsplit(infile, '_')[[1]][4:7], collapse = '_')
  outname = sprintf('abdb_figures_2019/%s_internet.pdf', name_parts)
  if (resize==T){
    outname = sprintf('abdb_figures_2019/%s_internet_resized.pdf', name_parts)
  }
  print(outname)
  # stop()
  df = read_csv(infile)#[1:100,]
  print(df)
  motifdf = tibble('source'=df$ab_motif, 'target' = df$ag_motif) %>% mutate('source_target' = paste0(source, '_', target))
  motifdf$target = sprintf('%s*', motifdf$target)
  print(length(unique(motifdf$source)))
  count_motifdf = count(motifdf, source, target, source_target)
  count_source = count(count_motifdf, source, sort = TRUE) %>% mutate('motif_source'= sprintf('Paratope', seq(1, nrow(.))))
  count_target = count(count_motifdf, target, sort = TRUE) %>% mutate('motif_source'= sprintf('Epitope', seq(1, nrow(.))))
  # print(sprintf('Paratope', seq(1,nrow(count_source))))
  print(count_source)
  print(count_target)
  print(count_target)
  print(motifdf)
  motif_source = c(sprintf('Paratope', seq(1, nrow(motifdf))), sprintf('Epitope', seq(1, nrow(motifdf))))
  source_target = data_frame('motif' = c(motifdf$source, motifdf$target), 'motif_source' = motif_source)
  print(source_target)
  st_motif_count = count(source_target, motif, motif_source, sort = TRUE)
  print(st_motif_count)
  # stop()
  nodes = data.frame(id = c(count_source$source, count_target$target), group = c(count_source$motif_source,count_target$motif_source),
                     score = as.integer(c(count_source$n, count_target$n)))
  nodes2 = data.frame(id=st_motif_count$motif, score = as.integer(st_motif_count$n), stringsAsFactors=FALSE, group=st_motif_count$motif_source)
  edges = data.frame(source = count_motifdf$source, target = count_motifdf$target, weight = count_motifdf$n, stringsAsFactors=FALSE)
  
  node_outfile = 'abdb_outfiles_2019/paratope_epitope_internet_nodes.csv'
  write.csv(x=nodes, file = node_outfile, row.names = FALSE)
  edge_outfile = 'abdb_outfiles_2019/paratope_epitope_internet_edges.csv'
  write.csv(x=edges, file = edge_outfile, row.names = FALSE)
  print('exiting...')
  # edges = data.frame(source = motifdf$source, target = motifdf$target, stringsAsFactors=FALSE)
  # ###
  # # Berfore bulding the network, plot degree distribution
  # ###
  # ### test for powerlaw
  # m_pl = displ$new(nodes$score)
  # est = estimate_xmin(m_pl)
  # print(est)
  # bs_p = bootstrap_p(m_pl) #boostraping for statistical significance according to Clauset et al. (2009)
  # print(bs_p)
  # power_label = sprintf('alpha: %s, p-value: %s', round(est$pars,1), round(bs_p$p,1))
  # ### end power law test
  # ### degree distribution
  # topn = 200
  # topndf = nodes[order(-nodes$score),][1:topn,]
  # topndf$group = factor(topndf$group, levels = c('Paratope', 'Epitope'))
  # ggplot(data = topndf[order(-topndf$score),]) +
  #   geom_bar(mapping = aes(x=reorder(id,-score), y=score, fill=group), stat = 'identity')+
  #   theme(axis.text.x  = element_text(angle = 90), 
  #         axis.title.y = element_text(size = 48),
  #         axis.title.x = element_text(size = 60),
  #         legend.text = element_text(size = 60),
  #         legend.title = element_blank(),
  #         panel.background = element_blank()) +
  #   # theme(axis.text.x  = element_text(angle = 90), axis.title = element_text(size = 28),
  #   #       legend.text = element_text(size = 28),
  #   #       legend.title = element_text(size = 28), 
  #   #       panel.background = element_blank())+
  #   labs(x='Interaction motifs (top 200 by degree)', y='# of edges (degree)', fill='Motif source')+
  #   scale_fill_manual(values = c(abcolor3, agcolor)) +
  #   geom_text(mapping = aes(x=id, y=score+4, label = score), angle = 90) + 
  #   geom_text(mapping = aes(x=180, y = 25, label = power_label), size = 10)
  # outpdf(outname, 'degree', width = 28, height = 10)
  # ### end degree distribution
  # ### cumulative degree distribution
  # print(topndf)
  # topndf = topndf %>% mutate(csum = cumsum(score)) %>% mutate(csum_percent = round(csum/sum(nodes$score),2))
  # csum_label = sprintf('Total number of edges: %s', sum(nodes$score))
  # ggplot(data = topndf) +
  #   geom_bar(mapping = aes(x=reorder(id,csum), y=csum, fill=group), stat = 'identity')+
  #   theme(axis.text.x  = element_text(angle = 90), 
  #         axis.title.y = element_text(size = 48),
  #         axis.title.x = element_text(size = 60),
  #         legend.text = element_text(size = 60),
  #         legend.title = element_blank()) +
  #   # theme(axis.text.x  = element_text(angle = 90), axis.title = element_text(size = 28),
  #   #       legend.text = element_text(size = 28),
  #   #       legend.title = element_text(size = 28)) +
  #         # panel.background = element_blank())+
  #   labs(x='Interaction motifs (top 200 by degree)', y='# of edges (degree) - percent', fill = 'Motif source')+
  #   scale_fill_manual(values = c(abcolor3, agcolor)) +
  #   geom_text(mapping = aes(x=id, y=csum+150, label = sprintf('%s - %s', csum, percent(csum_percent)), angle = 90)) + 
  #   geom_text(mapping = aes(x=20, y = 3300, label = csum_label), size = 10)
  # outpdf(outname, 'degree_cumsum', width = 28, height = 15)
  # print(topndf)
  # stop()
  # ### end cumulative degree distribution
  # ### end power law test
  # ### end the degree bit
  if(file.exists(outname)){
    response = file.remove(outname)
  }
  deleteAllNetworks()
  print(getNodeShapes())
  setNodeShapeDefault('ELLIPSE')
  createNetworkFromDataFrames(nodes, edges)
  # layoutNetwork('attribute-circle')
  print(getLayoutNames())
  node_attribute = getTableColumns()
  # print(node_attribute)
  node_size = as.numeric(node_attribute$score)
  node_attribute[node_attribute=='Paratope'] = '#3399CC'
  node_attribute[node_attribute=='Epitope'] = greys[5]
  node_attribute$group2 = node_attribute$group
  node_attribute$group2[node_attribute$group2==greys[5]] = greys[1]
  print(node_attribute[1000:2000,])
  # stop()
  if (resize==T){
    node_size = sqrt(as.numeric(node_attribute$score))*10
  }
  # stop()
  col_list = list("Paratope" = '#3399CC', "Epitope"= '#C0C0C0')
  print(head(node_attribute[order(-node_attribute$score),]))
  print(dim(node_attribute))
  print(head(nodes))
  print(dim(nodes))
  # stop()
  setNodeSizeBypass(node.names = node_attribute$name, new.sizes = node_size)
  setNodeColorBypass(node.names = node_attribute$name, new.colors = node_attribute$group)
  setNodeFillOpacityBypass(node.names = node_attribute$name, new.values = rep(200, nrow(node_attribute)))
  setNodeFontSizeBypass(node.names = node_attribute$name, new.sizes = node_size)
  setEdgeLineWidthDefault(new.width = 1)
  setEdgeColorDefault(new.color = greys[3])
  setNodeLabelColorBypass(node.names = node_attribute$name, new.colors = node_attribute$group)
  if (resize==T){
    setNodeColorBypass(node.names = node_attribute$name, new.colors = node_attribute$group2)
                            
  }
  bundleEdges()
  exportImage(outname, type = 'PDF')
  system(sprintf('open %s', outname))
  deleteAllNetworks()
}


paratope_epitope_motif_interaction_structures = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle.csv'
  name_parts = paste(strsplit(infile, '_')[[1]][4:7], collapse = '_')
  outname = sprintf('abdb_figures_2019/%s_internet_structure.pdf', name_parts)
  print(outname)
  df = read_csv(infile)#[1:100,]
  print(df)
  motifdf = tibble('source'=df$ab_motif, 'target' = df$ag_motif) %>% mutate('source_target' = paste0(source, '_', target))
  motifdf['pdbid'] = df$pdbid
  motifdf['segment'] = df$segment
  print(motifdf)
  motifdf$target = sprintf('%s*', motifdf$target)
  print(length(unique(motifdf$source)))
  print(motifdf)
  # stop()
  count_motifdf_structure = count(motifdf, source, target, source_target, pdbid, segment)
  count_motifdf = count(motifdf, source, target, source_target)
  count_source_structure = count(count_motifdf_structure, source, pdbid, segment, sort = FALSE) %>% mutate('motif_source'= sprintf('Paratope', seq(1, nrow(.))))
  count_source_structure = count(count_source_structure, source, pdbid, sort = F) %>% mutate('motif_source'= sprintf('Paratope', seq(1, nrow(.))))
  count_source_structure = count(count_source_structure, source, sort = F) %>% mutate('motif_source'= sprintf('Paratope', seq(1, nrow(.))))
  count_source = count(count_motifdf, source, sort = F) %>% mutate('motif_source'= sprintf('Paratope', seq(1, nrow(.))))
  count_target_structure = count(count_motifdf_structure, target, pdbid, segment, sort = F) %>% mutate('motif_source'= sprintf('Epitope', seq(1, nrow(.))))
  count_target_structure = count(count_target_structure, target, pdbid, sort = F) %>% mutate('motif_source'= sprintf('Epitope', seq(1, nrow(.))))
  count_target_structure = count(count_target_structure, target, sort = F) %>% mutate('motif_source'= sprintf('Epitope', seq(1, nrow(.))))
  count_target = count(count_motifdf, target, sort = F) %>% mutate('motif_source'= sprintf('Epitope', seq(1, nrow(.))))
  # print(sprintf('Paratope', seq(1,nrow(count_source))))
  print(count_source)
  print(count_source_structure)
  print(count_target)
  print(count_target_structure)
  # print(motifdf)
  ### replace degree with # of structures
  print(count_source)
  count_source$n = count_source_structure$n
  print(count_source)
  print(count_target)
  count_target$n = count_target_structure$n
  print(count_target)
  # stop()
  motif_source = c(sprintf('Paratope', seq(1, nrow(motifdf))), sprintf('Epitope', seq(1, nrow(motifdf))))
  source_target = data_frame('motif' = c(motifdf$source, motifdf$target), 'motif_source' = motif_source)
  print(source_target)
  st_motif_count = count(source_target, motif, motif_source, sort = TRUE)
  print(st_motif_count)
  # stop()
  nodes = data.frame(id = c(count_source$source, count_target$target), group = c(count_source$motif_source,count_target$motif_source),
                     score = as.integer(c(count_source$n, count_target$n)))
  nodes2 = data.frame(id=st_motif_count$motif, score = as.integer(st_motif_count$n), stringsAsFactors=FALSE, group=st_motif_count$motif_source)
  edges = data.frame(source = count_motifdf$source, target = count_motifdf$target, weight = count_motifdf$n, stringsAsFactors=FALSE)
  
  node_outfile = 'abdb_outfiles_2019/paratope_epitope_internet_nodes_structure.csv'
  write.csv(x=nodes, file = node_outfile, row.names = FALSE)
  edge_outfile = 'abdb_outfiles_2019/paratope_epitope_internet_edges_structure.csv'
  write.csv(x=edges, file = edge_outfile, row.names = FALSE)
  print('exiting...')
  # edges = data.frame(source = motifdf$source, target = motifdf$target, stringsAsFactors=FALSE)
  # ###
  # # Berfore bulding the network, plot degree distribution
  # ###
  # ### test for powerlaw
  # m_pl = displ$new(nodes$score)
  # est = estimate_xmin(m_pl)
  # print(est)
  # bs_p = bootstrap_p(m_pl) #boostraping for statistical significance according to Clauset et al. (2009)
  # print(bs_p)
  # power_label = sprintf('alpha: %s, p-value: %s', round(est$pars,1), round(bs_p$p,1))
  # ### end power law test
  # ### degree distribution
  # topn = 200
  # topndf = nodes[order(-nodes$score),][1:topn,]
  # topndf$group = factor(topndf$group, levels = c('Paratope', 'Epitope'))
  # ggplot(data = topndf[order(-topndf$score),]) +
  #   geom_bar(mapping = aes(x=reorder(id,-score), y=score, fill=group), stat = 'identity')+
  #   theme(axis.text.x  = element_text(angle = 90), axis.title = element_text(size = 28),
  #         legend.text = element_text(size = 28),
  #         legend.title = element_text(size = 28),
  #         panel.background = element_blank())+
  #   labs(x='Interaction motifs (top 200 by degree)', y='# of edges (degree)', fill='Motif source')+
  #   scale_fill_manual(values = c(abcolor3, agcolor)) +
  #   geom_text(mapping = aes(x=id, y=score+4, label = score), angle = 90) +
  #   geom_text(mapping = aes(x=180, y = 25, label = power_label), size = 10)
  # outpdf(outname, 'degree', width = 28, height = 10)
  # ### end degree distribution
  # ### cumulative degree distribution
  # print(topndf)
  # topndf = topndf %>% mutate(csum = cumsum(score)) %>% mutate(csum_percent = round(csum/sum(nodes$score),2))
  # csum_label = sprintf('Total number of edges: %s', sum(nodes$score))
  # ggplot(data = topndf) +
  #   geom_bar(mapping = aes(x=reorder(id,csum), y=csum, fill=group), stat = 'identity')+
  #   theme(axis.text.x  = element_text(angle = 90), axis.title = element_text(size = 28),
  #         legend.text = element_text(size = 28),
  #         legend.title = element_text(size = 28)) +
  #         # panel.background = element_blank())+
  #   labs(x='Interaction motifs (top 200 by degree)', y='# of edges (degree) - percent', fill = 'Motif source')+
  #   scale_fill_manual(values = c(abcolor3, agcolor)) +
  #   geom_text(mapping = aes(x=id, y=csum+150, label = sprintf('%s - %s', csum, percent(csum_percent)), angle = 90)) +
  #   geom_text(mapping = aes(x=20, y = 3300, label = csum_label), size = 10)
  # outpdf(outname, 'degree_cumsum', width = 28, height = 15)
  # print(topndf)
  # # stop()
  # ### end cumulative degree distribution
  # ### end power law test
  # ### end the degree bit
  if(file.exists(outname)){
    response = file.remove(outname)
  }
  deleteAllNetworks()
  print(getNodeShapes())
  setNodeShapeDefault('ELLIPSE')
  createNetworkFromDataFrames(nodes, edges)
  # layoutNetwork('attribute-circle')
  print(getLayoutNames())
  node_attribute = getTableColumns()
  node_attribute = node_attribute[order(node_attribute$name),]
  print(head(node_attribute))
  print(head(nodes[order(nodes$id),]))
  node_size = sqrt(as.numeric(node_attribute$score))*7
  # print(node_size)
  # stop()
  col_list = list("Paratope" = '#3399CC', "Epitope"= '#C0C0C0')
  node_attribute[node_attribute=='Paratope'] = '#3399CC'
  node_attribute[node_attribute=='Epitope'] = greys[5]
  print(head(node_attribute[order(-node_attribute$score),]))
  print(dim(node_attribute))
  print(head(nodes))
  print(dim(nodes))
  # stop()
  setNodeSizeBypass(node.names = node_attribute$name, new.sizes = node_size)
  setNodeColorBypass(node.names = node_attribute$name, new.colors = node_attribute$group)
  setNodeFillOpacityBypass(node.names = node_attribute$name, new.values = rep(200, nrow(node_attribute)))
  setNodeFontSizeBypass(node.names = node_attribute$name, new.sizes = node_size)
  setEdgeLineWidthDefault(new.width = 1)
  setEdgeColorDefault(new.color = greys[3])
  setNodeLabelColorBypass(node.names = node_attribute$name, new.colors = node_attribute$group)
  bundleEdges()
  exportImage(outname, type = 'PDF')
  system(sprintf('open %s', outname))
  deleteAllNetworks()
}



random_internet = function(){
  infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle.csv'
  df = read_csv(infile)
  df$ag_motif = sprintf('%s*', df$ag_motif)
  nsamples = 1000
  nrepeats = 100
  dfs = c(seq(1:nrepeats))
  # print(unique(df$ab_motif))
  # stop()
  abmotiflist = lapply(unique(df$ab_motif), function(x) x=0)
  names(abmotiflist)  = unique(df$ab_motif)
  agmotiflist = lapply(unique(df$ag_motif), function(x) x=0)
  names(agmotiflist)  = unique(df$ag_motif)
  reps = paste0('rep', seq(1:nrepeats))
  stdf = data_frame()
  print(stdf)
  for (i in dfs){
    print(i)
    abmotifs= df$ab_motif %>% sample(., nsamples, replace = TRUE)
    agmotifs= df$ag_motif %>% sample(., nsamples, replace = TRUE)
    print(agmotifs)
    # rabmotifs = sample(abmotifs, 100)
    rdf = data_frame('source' = abmotifs, 'target' = agmotifs, 'rep' = rep(i))
    count_rdf = count(rdf, source, target, rep, sort = TRUE)
    print(count_rdf)
    # rbind(sourcelist[reps[i]],count_rdf)
    stdf = rbind(stdf, count_rdf)
  }
  edges = stdf %>% group_by(source, target) %>% summarise(n = median(n))
  source_nodes = stdf %>% count(., source, rep, sort = TRUE) %>% group_by(source) %>% 
    summarise(n = round(mean(n))) %>% mutate('group' = rep('Paratope', nrow(.))) %>%
    mutate(n_reps = stdf %>% count(source, rep) %>% group_by(source) %>% summarise(l = length(source)) %>% .$l ) %>%
    mutate(n_samples = stdf %>% group_by(source) %>% summarise(l = length(source)) %>% .$l ) %>%
    mutate(sd = stdf %>% count(source, rep) %>% group_by(source) %>% summarise(sd = sd(n)) %>% .$sd ) %>%
    mutate(ste_mean = sd/sqrt(n_reps))
  target_nodes = stdf %>% count(., target, rep, sort = TRUE) %>% group_by(target) %>% 
    summarise(n = round(median(n))) %>% mutate('group' = rep('Epitope', nrow(.))) %>%
    mutate(n_reps = stdf %>% count(target, rep) %>% group_by(target) %>% summarise(l = length(target)) %>% .$l ) %>%
    mutate(n_samples = stdf %>% group_by(target) %>% summarise(l = length(target)) %>% .$l ) %>%
    mutate(sd = stdf %>% count(target, rep) %>% group_by(target) %>% summarise(sd = sd(n)) %>% .$sd ) %>%
    mutate(ste_mean = sd/sqrt(n_reps))
  print(target_nodes[order(-target_nodes$n),])
  print(source_nodes[order(-source_nodes$n),])
  colnames(source_nodes)[1] = 'motif'
  colnames(target_nodes)[1] = 'motif'
  print(source_nodes)
  print(target_nodes)
  outdf = rbind(source_nodes, target_nodes) # remember my child this is the nodes
  ### test for powerlaw
  outdf_nozeros = outdf[outdf$n >0,]
  print(outdf_nozeros)
  m_pl = displ$new(outdf_nozeros$n)
  est = estimate_xmin(m_pl)
  print(est)
  bs_p = bootstrap_p(m_pl) #boostraping for statistical significance according to Clauset et al. (2009)
  print(bs_p)
  power_label = sprintf('alpha: %s, p-value: %s', round(est$pars,1), round(bs_p$p,1))
  print(power_label)
  ### end power law test
  # stop()
  topn = 200
  topndf = outdf[order(-outdf$n),][1:topn,]
  topndf[topndf=='Paratope'] = 'Randomly sampled paratope'
  topndf[topndf=='Epitope'] = 'Randomly sampled epitope'
  topndf$group = factor(topndf$group, levels = c('Randomly sampled paratope', 'Randomly sampled epitope'))
  # stop()
  ### degree distribution
  ggplot(data = topndf) +
    geom_bar(mapping = aes(x=reorder(motif, -n), y=n, fill = group), stat = 'identity')+
    geom_text(mapping = aes(x=motif, y=n+ste_mean+4, label = n), angle=90) +
    theme(axis.text.x  = element_text(angle = 90),
          axis.title.y = element_text(size = 48),
          axis.title.x = element_text(size = 60),
          legend.text = element_text(size = 60),
          legend.title = element_blank(),
          panel.background = element_blank()) +
    # theme(axis.text.x  = element_text(angle = 90), axis.title = element_text(size = 28),
    #       legend.text = element_text(size = 28),
    #       legend.title = element_text(size = 28),
          # panel.background = element_blank())+
    labs(x='Interaction motifs (top 200 by degree)', y='# of edges \n(mean degree)', fill='Motif source') +
    scale_fill_manual(values = c(rabcolor, ragcolor)) +
    geom_text(mapping = aes(x=186, y=15, label=power_label), size=10)+
    geom_errorbar(mapping = aes(x= motif, y=n, ymin=n-ste_mean, ymax= n+ ste_mean), size=0.5, color=my_spectral[2])
  outpdf(infile, 'random_internet_degree', width = 28, height = 10 )
  ### end degree distribution
  ### cumulative degree distribution
  topndf = topndf %>% mutate(csum = cumsum(n)) %>% mutate(csum_percent = round(csum/sum(outdf$n),2))
  print(topndf)
  csum_label = sprintf('Total number of edges: %s', sum(outdf$n))
  print(tail(topndf))
  print(csum_label)
  ggplot(data = topndf) +
    geom_bar(mapping = aes(x=reorder(motif,csum), y=csum, fill=group), stat = 'identity')+
    theme(axis.text.x  = element_text(angle = 90),
          axis.title.y = element_text(size = 48),
          axis.title.x = element_text(size = 60),
          legend.text = element_text(size = 60),
          legend.title = element_blank()) +
    # theme(axis.text.x  = element_text(angle = 90), axis.title = element_text(size = 28),
    #       legend.text = element_text(size = 28),
    #       legend.title = element_text(size = 28))+
    labs(x='Interaction motifs (top 200 by degree)', y='# of edges (mean degree) - percent', fill = 'Motif source')+
    scale_fill_manual(values = c(rabcolor, ragcolor)) +
    geom_text(mapping = aes(x=motif, y=csum+150, label = sprintf('%s - %s', csum, percent(csum_percent)), angle = 90)) +
    geom_text(mapping = aes(x=20, y = 3300, label = csum_label), size = 10)
  outpdf(infile, 'random_internet_degree_cumsum', width = 28, height = 15)
  ### end cum degree distribution
  stop()
  print(topndf)
  ###
  # network bit
  ###
  name_parts = paste(strsplit(infile, '_')[[1]][4:7], collapse = '_')
  outname = sprintf('abdb_figures_2019/%s_random_internet.pdf', name_parts)
  motifdf = stdf %>% mutate(source_target = sprintf(stdf$source, '_', stdf$target))
  count_motifdf = count(motifdf, source, target, source_target)
  count_source = count(count_motifdf, source, sort = TRUE) %>% mutate('motif_source'= sprintf('Paratope', seq(1, nrow(.))))
  count_target = count(count_motifdf, target, sort = TRUE) %>% mutate('motif_source'= sprintf('Epitope', seq(1, nrow(.))))
  print(count_motifdf)
  motif_source = c(sprintf('Paratope', seq(1, nrow(motifdf))), sprintf('Epitope', seq(1, nrow(motifdf))))
  source_target = data_frame('motif' = c(motifdf$source, motifdf$target), 'motif_source' = motif_source)
  # print(source_target)
  st_motif_count = count(source_target, motif, motif_source, sort = TRUE)
  # print(st_motif_count)
  print(count_motifdf)
  # stop()
  # nodes = data.frame(id = c(count_source$source, count_target$target), group = c(count_source$motif_source,count_target$motif_source),
  #                    score = as.integer(c(count_source$n, count_target$n)))
  nodes = data.frame(id = outdf$motif, group = outdf$group, score = outdf$n)
  # edges = data.frame(source = count_motifdf$source, target = count_motifdf$target, weight = count_motifdf$n, stringsAsFactors=FALSE)
  edges = data.frame(source = edges$source, target = edges$target, weight = edges$n, stringsAsFactors=FALSE)
  # edges = edges[1:10,]
  print(head(edges))
  print(head(nodes))
  node_outfile = 'abdb_outfiles_2019/random_internet_nodes.csv'
  write.csv(x=nodes, file = node_outfile, row.names = FALSE)
  edge_outfile = 'abdb_outfiles_2019/random_internet_edges.csv'
  write.csv(x=edges, file = edge_outfile, row.names = FALSE)
  print('exiting...')
  # stop()
  if(file.exists(outname)){
    response = file.remove(outname)
  }
  deleteAllNetworks()
  print(getNodeShapes())
  setNodeShapeDefault('ELLIPSE')
  createNetworkFromDataFrames(nodes, edges)
  # layoutNetwork('attribute-circle')
  print(getLayoutNames())
  node_attribute = getTableColumns()
  # print(node_attribute)
  node_size = as.numeric(node_attribute$score)
  node_size = sqrt(node_size)*5
  # stop()
  col_list = list("Paratope" = '#3399CC', "Epitope"= '#C0C0C0')
  node_attribute[node_attribute=='Paratope'] = rabcolor
  node_attribute[node_attribute=='Epitope'] = ragcolor
  print(head(node_attribute[order(-node_attribute$score),]))
  print(dim(node_attribute))
  print(head(nodes))
  print(dim(nodes))
  # stop()
  setNodeSizeBypass(node.names = node_attribute$name, new.sizes = node_size)
  setNodeColorBypass(node.names = node_attribute$name, new.colors = node_attribute$group)
  setNodeFillOpacityBypass(node.names = node_attribute$name, new.values = rep(200, nrow(node_attribute)))
  setNodeFontSizeBypass(node.names = node_attribute$name, new.sizes = node_size)
  setEdgeLineWidthDefault(new.width = 1)
  setEdgeColorDefault(new.color = greys[3])
  setNodeLabelColorBypass(node.names = node_attribute$name, new.colors = node_attribute$group)
  bundleEdges()
  exportImage(outname, type = 'PDF')
  system(sprintf('open %s', outname))
  deleteAllNetworks()
  ###
}

degree_correlation = function(){
  rinfile = 'abdb_outfiles_2019/random_internet_nodes.csv'
  peinfile = 'abdb_outfiles_2019/paratope_epitope_internet_nodes.csv'
  rdf = read_csv(rinfile)
  cutoff = 10
  rdf = rdf[rdf$score >=cutoff,]
  inmotifs = rdf$id 
  print(rdf)
  pedf =read_csv(peinfile)
  pedf = pedf[pedf$id %in% inmotifs,]
  print(rdf)
  print(pedf)
  mdf = merge(rdf, pedf, by='id')
  print(head(mdf))
  topnmdf = mdf[order(-mdf$score.x),][1:nrow(mdf),]
  scor = cor(topnmdf$score.x,topnmdf$score.y, method = 'spearman')
  scor_label = sprintf('Spearman cor.: %s\nn: %s', round(scor,2), nrow(mdf))
  print(scor_label)
  # stop()
  ggplot(data=topnmdf) + 
    geom_point(mapping = aes(x=score.x, y=score.y), color = alpha(my_spectral[3],0.4), size=5) +
    geom_smooth(mapping = aes(x=score.x, y=score.y), color = alpha(my_spectral[3], 0.9), method = 'lm') + 
    geom_text(mapping = aes(x=0, y=170, label = scor_label), size=10, hjust=0) + 
    geom_density_2d(mapping = aes(x=score.x, y=score.y), color = my_spectral[3]) +
    labs(x=sprintf('Randomly sampled \n(degree >= %s)', cutoff), y='Observed node degree') + 
    theme(axis.title = element_text(size=20), panel.background = element_blank()) + 
    theme(axis.title = element_text(size=35),
          axis.text = element_text(size=20),
          panel.background = element_blank())
  outpdf('random_paraepi_degree', 'corr')
}


cross_reactivity = function() {
  infile = 'abdb_outfiles_2019/paratope_epitope_internet_edges_paratope_cross_top10.csv'
  df = read_csv(infile)
  # df = df[1:100, 1:100]
  print(df)
  dfmat = df[,-1]
  print(dfmat)
  name_parts = strsplit(infile, '\\.')[[1]][1]
  outname = sprintf('%s_heatmap.pdf', name_parts)
  print(outname)
  # stop()
  hmap = pheatmap(dfmat,
                  labels_row = df$motif,
                  # color = my_spectral[1:14],
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  color = viridis(7),
                  number_color = 'gray',
                  cellwidth = 11, cellheight = 11,
                  border_color = NA,
                  # display_numbers = TRUE,
                  # cutree_rows = 4,
                  legend = TRUE,
                  # legend_breaks = c(0.2,0.4,0.6,0.8,1),
                  # legend_labels =  c('0.2', '0.4', '0.6', '0.8','Sequential frequency\n1')
                  )
  pdf(outname, width = 8, height = 8)
  grid.draw(hmap$gtable)
  dev.off()
  system(sprintf('open %s', outname))
  
}


cross_reactivity_density = function() {
  infile = 'abdb_outfiles_2019/paratope_epitope_internet_edges_paratope_cross_density.csv'
  df = read_csv(infile)
  pdf = subset(df, motif_source=='paratope')
  rpdf = subset(df, motif_source =='rparatope')
  edf = subset(df, motif_source=='epitope')
  redf = subset(df, motif_source =='repitope')
  print(pdf)
  print(rpdf)
  pks = ks.test(pdf$percent_overlap, rpdf$percent_overlap)
  eks = ks.test(edf$percent_overlap, redf$percent_overlap)
  print(str(pks))
  print(pks$statistic)
  pks_label = sprintf('KS test; D: %s; p-value: %s', round(pks$statistic,2), '< 10e-10')
  eks_label = sprintf('KS test; D: %s; p-value: %s', round(eks$statistic,2), '< 10e-10')
  print(pks_label)
  print(eks_label)
  ksdf = tibble('motif_source' = c('paratope', 'epitope'), 'kstest' = c(pks_label, eks_label))
  print(ksdf)
  df[df=='paratope'] = 'Paratope'
  df[df=='epitope'] = 'Epitope'
  df[df=='repitope'] = 'Randomly sampled epitope'
  df[df=='rparatope'] = 'Randomly sampled paratope'
  df$motif_source  = factor(df$motif_source, levels = c('Paratope', 'Epitope', 'Randomly sampled paratope', 'Randomly sampled epitope'))
  meds  = df %>% group_by(motif_source) %>% summarise(med = median(percent_overlap))
  print(df)
  med_labels = sprintf('Median: %s', percent(meds$med/100))
  print(meds)
  #Function to transform data to y positions
  break_point = 0.05
  trans <- function(x){pmin(x,break_point) + 0.05*pmax(x-break_point,0)}
  yticks= c(0.01,0.02, 0.03, 0.04, 0.05, 0.2, 0.3, 0.4, 0.5)
  yticks_labels = c('0.01','0.02', '0.03', '0.04', '\\ \\', '0.2', '0.3', '0.4', '0.5')
  print(trans(yticks))
  # stop()
  ggplot(data= df) + 
    geom_density(mapping = aes(x=percent_overlap, y=trans(..density..), fill=motif_source, color = motif_source), color = NA, alpha = 0.6) + 
    geom_text(data = meds, mapping = aes(x=c(2, 2, 14, 14), y = trans(c(0.5, 0.4, 0.5, 0.4)), label= med_labels, color=motif_source), hjust=0, show.legend = FALSE) + 
    geom_text(data = ksdf, mapping = aes(x=c(26, 26), y = trans(c(0.5, 0.4)), label= c(pks_label, eks_label)), hjust =0) + 
    scale_fill_manual(values = c(abcolor3, agcolor, rabcolor, ragcolor)) +  
    scale_color_manual(values = c(abcolor3, agcolor, rabcolor, ragcolor)) + 
    scale_y_continuous(limits=c(0,NA), breaks=trans(yticks), labels=yticks_labels) +
    geom_rect(aes(xmin=0, xmax=2, ymin=trans(0.05), ymax=trans(0.06)), fill="white") + 
    labs(x= 'Interaction partner overlap (%)', y='Density', fill = 'Motif source', color = 'Motif source') + 
    theme(legend.position = 0,
          axis.title = element_text(size = 30),
          axis.text = element_text(size = 20),
          strip.text = element_text(size = 12))
  outpng('motif_partner', 'density', width = 20)

}


cross_reactivity_density_logscale = function() {
  infile = 'abdb_outfiles_2019/paratope_epitope_internet_edges_paratope_cross_density.csv'
  df = read_csv(infile)
  pdf = subset(df, motif_source=='paratope')
  rpdf = subset(df, motif_source =='rparatope')
  edf = subset(df, motif_source=='epitope')
  redf = subset(df, motif_source =='repitope')
  print(pdf)
  print(rpdf)
  pks = ks.test(pdf$percent_overlap, rpdf$percent_overlap)
  eks = ks.test(edf$percent_overlap, redf$percent_overlap)
  print(str(pks))
  print(pks$statistic)
  pks_label = sprintf('KS test; D: %s; p-value: %s', round(pks$statistic,2), '< 10e-10')
  eks_label = sprintf('KS test; D: %s; p-value: %s', round(eks$statistic,2), '< 10e-10')
  print(pks_label)
  print(eks_label)
  ksdf = tibble('motif_source' = c('paratope', 'epitope'), 'kstest' = c(pks_label, eks_label))
  print(ksdf)
  df[df=='paratope'] = 'Paratope'
  df[df=='epitope'] = 'Epitope'
  df[df=='repitope'] = 'Random epitope'
  df[df=='rparatope'] = 'Random paratope'
  df$motif_source  = factor(df$motif_source, levels = c('Paratope', 'Epitope', 'Random paratope', 'Random epitope'))
  meds  = df %>% group_by(motif_source) %>% summarise(med = median(percent_overlap))
  print(df)
  med_labels = sprintf('Median: %s', percent(meds$med/100))
  print(meds)
  ggplot(data= df) + 
    geom_density(mapping = aes(x=percent_overlap, y=-1/log(..density..), fill=motif_source, color = motif_source), color = NA, alpha = 0.6) + 
    geom_text(data = meds, mapping = aes(x=c(5, 5, 16, 16), y = c(1.2, 1, 1.2, 1), label= med_labels, color=motif_source), hjust=0) + 
    geom_text(data = ksdf, mapping = aes(x=c(28, 28), y = c(1.2, 1), label= c(pks_label, eks_label)), hjust =0) + 
    scale_fill_manual(values = c(abcolor3, agcolor, rabcolor, ragcolor)) +  
    scale_color_manual(values = c(abcolor3, agcolor, rabcolor, ragcolor)) + 
    labs(x= 'Interaction partner overlap (%)', y='Density (-1/log_density)', fill = 'Motif source', color = 'Motif source')
  outpdf('motif_partner', 'density_log', width = 20)

}

cross_reactivity_degree_cor = function(){
  rinfile = 'abdb_outfiles_2019/random_internet_nodes.csv'
  peinfile = 'abdb_outfiles_2019/paratope_epitope_internet_nodes.csv'
  rdf = read_csv(rinfile)
  rdf[rdf == 'Paratope'] = 'rparatope'
  rdf[rdf == 'Epitope'] = 'repitope'
  pedf = read_csv(peinfile)
  pedf[pedf == 'Paratope'] = 'rparatope'
  pedf[pedf == 'Epitope'] = 'repitope'
  print(rdf)
  print(pedf)
  mdegre_df = rbind(pedf, rdf)
  print(mdegre_df)
  infile = 'abdb_outfiles_2019/paratope_epitope_internet_edges_paratope_cross_density.csv'
  df = read_csv(infile) %>%.[.$percent_overlap < 100,] # excludes self pairing
  meancross = df %>% group_by(.,motif1, motif_source) %>% summarise('ave' = mean(percent_overlap)) %>% .[order(-.$ave),]
  colnames(meancross)[1] = 'id'
  merge_degree_overlapdf = merge(mdegre_df,meancross, by = 'id')
  print(meancross)
  print(head(merge_degree_overlapdf))
  lvls  = c('paratope', 'epitope', 'rparatope', 'repitope')
  facet_labels = c('Paratope', 'Epitope', 'Randomly sampled paratope', 'Randomly sampled epitope')
  names(facet_labels) = lvls
  merge_degree_overlapdf$motif_source = factor(merge_degree_overlapdf$motif_source, levels = lvls)
  for (msource in unique(merge_degree_overlapdf$motif_source)) {merge_degree_overlapdf[]}
  ggplot(data =merge_degree_overlapdf) + 
    geom_point(mapping = aes(x=score, y=ave, color=motif_source), size = 5, alpha=0.5) + 
    facet_wrap(~ motif_source, ncol = 4, labeller = labeller(motif_source=facet_labels), scales = 'free') +
    scale_color_manual(values = c(abcolor3, agcolor, rabcolor, ragcolor)) + 
    labs(x='Node degree', y = 'Mean interaction \npartner overlap (%)') + 
    theme(axis.title = element_text(size = 25),
          legend.text = element_text(size = 25),
          axis.text = element_text(size = 25),
          legend.title = element_blank(),
          strip.text = element_text(size=15))
  outpdf('degree', 'mean_overlap_scatter', width = 15, height = 5)
  
  
}
  
# run stuff
# note to self: cytoscape style needs to be reset from rcy3 (can't find the switch for this) or restart the program
# note to self: must use dataframe, tibble would not work wit cytoscape
# paratope_epitope_motif_interaction()
# random_internet()
# degree_correlation()
# cross_reactivity()
# cross_reactivity_density()
# cross_reactivity_density_logscale()
# cross_reactivity_degree_cor()
# paratope_epitope_motif_interaction_structures()
paratope_epitope_motif_interaction(resize = T)



