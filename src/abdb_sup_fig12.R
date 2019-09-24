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

ppi_mcolor_rgb = c(col2rgb('orange'))/255
ppi_mcolor = rgb(ppi_mcolor_rgb[1], ppi_mcolor_rgb[2], ppi_mcolor_rgb[3])
ppi_pcolor_rgb = c(col2rgb('darkorange'))/255
ppi_pcolor = rgb(ppi_pcolor_rgb[1], ppi_pcolor_rgb[2], ppi_pcolor_rgb[3])


rppi_mcolor_rgb = c(col2rgb('darkolivegreen'))/255
rppi_mcolor = rgb(rppi_mcolor_rgb[1], ppi_mcolor_rgb[2], ppi_mcolor_rgb[3])
rppi_pcolor_rgb = c(col2rgb('darkolivegreen3'))/255
rppi_pcolor = rgb(rppi_pcolor_rgb[1], ppi_pcolor_rgb[2], ppi_pcolor_rgb[3])

random_mcolor = '#3399CC'
random_pcolor = '#C0C0C0'
print(sprintf('paratope: %s', abcolor3))
print(sprintf('epitope: %s', agcolor))
print(sprintf('ppi motif: %s', ppi_mcolor))
print(sprintf('ppi partner: %s', ppi_pcolor))
print(sprintf('radom motif: %s', random_mcolor))
print(sprintf('radom partner: %s', random_pcolor))

stop()

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



ppi_motif_interaction = function(){
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv'
  name_parts = paste(strsplit(infile, '_')[[1]][4:7], collapse = '_')
  outname = sprintf('abdb_figures_2019/%s_internet.pdf', name_parts)
  print(outname)
  # df = read_csv(infile)[1:100,]
  df = read_csv(infile)
  print(df)
  df1 = df[df$intertype1=='interdomain',]
  df2 = df[df$intertype1=='intradomain',]
  print(dim(df1))
  print(dim(df2))
  motifdf = tibble('source'=df1$gap_pattern1, 'target' = df1$gap_pattern2) %>% mutate('source_target' = paste0(source, '_', target))
  motifdf$target = sprintf('%s*', motifdf$target)
  print(length(unique(motifdf$source)))
  count_motifdf = count(motifdf, source, target, source_target)
  count_source = count(count_motifdf, source, sort = TRUE) %>% mutate('motif_source'= sprintf('PPI_motif', seq(1, nrow(.))))
  count_target = count(count_motifdf, target, sort = TRUE) %>% mutate('motif_source'= sprintf('PPI_partner', seq(1, nrow(.))))
  # print(sprintf('Paratope', seq(1,nrow(count_source))))
  print(count_source)
  print(count_target)
  print(motifdf)
  motif_source = c(sprintf('PPI_motif', seq(1, nrow(motifdf))), sprintf('PPI_partner', seq(1, nrow(motifdf))))
  source_target = data_frame('motif' = c(motifdf$source, motifdf$target), 'motif_source' = motif_source)
  print(tail(source_target))
  st_motif_count = count(source_target, motif, motif_source, sort = TRUE)
  print(st_motif_count)
  nodes = data.frame(id = c(count_source$source, count_target$target), group = c(count_source$motif_source,count_target$motif_source),
                     score = as.integer(c(count_source$n, count_target$n)))
  # nodes2 = data.frame(id=st_motif_count$motif, score = as.integer(st_motif_count$n), stringsAsFactors=FALSE, group=st_motif_count$motif_source)
  edges = data.frame(source = count_motifdf$source, target = count_motifdf$target, weight = count_motifdf$n, stringsAsFactors=FALSE)
  node_outfile = 'abdb_outfiles_2019/ppi_internet_nodes.csv'
  write.csv(x=nodes, file = node_outfile, row.names = FALSE)
  edge_outfile = 'abdb_outfiles_2019/ppi_internet_edges.csv'
  write.csv(x=edges, file = edge_outfile, row.names = FALSE)
  ###
  # Berfore bulding the network, plot degree distribution
  ###
  ### test for powerlaw
  m_pl = displ$new(nodes$score)
  est = estimate_xmin(m_pl)
  print(est)
  bs_p = bootstrap_p(m_pl) #boostraping for statistical significance according to Clauset et al. (2009)
  bs_p_label = round(bs_p$p,2)
  if (round(bs_p_label,10) == 0){bs_p_label = '<10e-10'}
  print(bs_p_label)
  power_label = sprintf('alpha: %s, p-value (power law): %s', round(est$pars,1), bs_p_label)
  ### end power law test
  ### degree distribution
  topn = 200
  topndf = nodes[order(-nodes$score),][1:topn,]
  # topndf$group = factor(topndf$group, levels = c('Paratope', 'Epitope'))
  ggplot(data = topndf[order(-topndf$score),]) +
    geom_bar(mapping = aes(x=reorder(id,-score), y=score, fill=group), stat = 'identity')+
    theme(axis.text.x  = element_text(angle = 90), axis.title = element_text(size = 28),
          legend.text = element_text(size = 28),
          legend.title = element_text(size = 28), 
          panel.background = element_blank())+
    labs(x='Interaction motifs (top 200 by degree)', y='# of edges (degree)', fill='Motif source')+
    scale_fill_manual(values = c(ppi_mcolor, ppi_pcolor)) +
    geom_text(mapping = aes(x=id, y=score+7, label = score), angle = 90) + 
    geom_text(mapping = aes(x=170, y = 50, label = power_label), size = 10)
  outpdf(outname, 'ppi_degree', width = 28, height = 10)
  ### end degree distribution
  ### cumulative degree distribution
  print(topndf)
  topndf = topndf %>% mutate(csum = cumsum(score)) %>% mutate(csum_percent = round(csum/sum(nodes$score),2))
  csum_label = sprintf('Total number of edges: %s', sum(nodes$score))
  csum_label_node = sprintf('Total number of nodes: %s', nrow(nodes))
  print(csum_label)
  print(csum_label_node)
  ggplot(data = topndf) +
    geom_bar(mapping = aes(x=reorder(id,csum), y=csum, fill=group), stat = 'identity')+
    theme(axis.text.x  = element_text(angle = 90), axis.title = element_text(size = 28),
          legend.text = element_text(size = 28),
          legend.title = element_text(size = 28)) +
          # panel.background = element_blank())+
    labs(x='Interaction motifs (top 200 by degree)', y='# of edges (degree) - percent', fill = 'Motif source')+
    scale_fill_manual(values = c(ppi_mcolor, ppi_pcolor)) +
    geom_text(mapping = aes(x=id, y=csum+500, label = sprintf('%s - %s', csum, percent(csum_percent)), angle = 90)) + 
    geom_text(mapping = aes(x=20, y = 8000, label = csum_label), size = 10) +
    geom_text(mapping = aes(x=20, y = 8000-500, label = csum_label_node), size = 10)  
  outpdf(outname, 'ppi_degree_cumsum', width = 28, height = 15)
  print(topndf)
  # stop()
  ### end cumulative degree distribution
  ### end power law test
  ### end the degree bit
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
  # stop()
  node_attribute[node_attribute=='PPI_motif'] = ppi_mcolor
  node_attribute[node_attribute=='PPI_partner'] = ppi_pcolor
  print(head(node_attribute[order(-node_attribute$score),]))
  print(dim(node_attribute))
  print(head(nodes))
  print(dim(nodes))
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



downsampled_ppi_motif_internet = function(){
  infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv'
  df = read_csv(infile)
  print(df)
  #in ppi case, abmotif is gap_patern1, and agmotif is gap_pattern2
  df$gap_pattern2 = sprintf('%s*', df$gap_pattern2)
  nsamples = 5000
  nrepeats = 100
  dfs = c(seq(1:nrepeats))
  # print(unique(df$ab_motif))
  # stop()
  #in ppi case, abmotif is gap_patern1, and agmotif is gap_pattern2
  abmotiflist = lapply(unique(df$gap_pattern1), function(x) x=0)
  names(abmotiflist)  = unique(df$gap_pattern1)
  agmotiflist = lapply(unique(df$gap_pattern2), function(x) x=0)
  names(agmotiflist)  = unique(df$gap_pattern2)
  reps = paste0('rep', seq(1:nrepeats))
  stdf = data_frame()
  print(stdf)
  for (i in dfs){
    print(i)
    # abmotifs= df$ab_motif %>% sample(., nsamples, replace = TRUE)
    # agmotifs= df$ag_motif %>% sample(., nsamples, replace = TRUE)
    dsdf = sample_n(df, nsamples, replace = TRUE)
    # rabmotifs = sample(abmotifs, 100)
    rdf = data_frame('source' = dsdf$gap_pattern1, 'target' = dsdf$gap_pattern2, 'rep' = rep(i))
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
  bs_p_label = round(bs_p$p,2)
  if (round(bs_p_label,10) == 0){bs_p_label = '<10e-10'}
  print(bs_p_label)
  power_label = sprintf('alpha: %s, p-value (power law): %s', round(est$pars,1), bs_p_label)
  print(power_label)
  ### end power law test
  # stop()
  topn = 200
  topndf = outdf[order(-outdf$n),][1:topn,]
  topndf[topndf=='Paratope'] = 'Downsampled PPI motif'
  topndf[topndf=='Epitope'] = 'Downsampled PPI motif partner'
  # topndf$group = factor(topndf$group, levels = c('Randomly sampled motif', 'Randomly sampled motif partner'))
  print(topndf)
  # stop()
  ### degree distribution
  ggplot(data = topndf) +
    geom_bar(mapping = aes(x=reorder(motif, -n), y=n, fill = group), stat = 'identity')+
    geom_text(mapping = aes(x=motif, y=n+ste_mean+4, label = n), angle=90) +
    theme(axis.text.x  = element_text(angle = 90), axis.title = element_text(size = 28),
          legend.text = element_text(size = 28),
          legend.title = element_text(size = 28),
          panel.background = element_blank()) +
    labs(x='Interaction motifs (top 200 by degree)', y='# of edges (mean degree)', fill='Motif source') +
    scale_fill_manual(values = c(rppi_mcolor, rppi_pcolor)) +
    geom_text(mapping = aes(x=170, y=18, label=power_label), size=10)+
    geom_errorbar(mapping = aes(x= motif, y=n, ymin=n-ste_mean, ymax= n+ ste_mean), size=0.5, color=my_spectral[2])
  outpdf(infile, 'downsampled_ppi_internet_degree', width = 28, height = 10 )
  ### end degree distribution
  ### cumulative degree distribution
  topndf = topndf %>% mutate(csum = cumsum(n)) %>% mutate(csum_percent = round(csum/sum(outdf$n),2))
  print(topndf)
  csum_label = sprintf('Total number of edges: %s', sum(outdf$n))
  csum_label_node = sprintf('Total number of nodes: %s', nrow(outdf))
  print(tail(topndf))
  print(csum_label)
  ggplot(data = topndf) +
    geom_bar(mapping = aes(x=reorder(motif,csum), y=csum, fill=group), stat = 'identity')+
    theme(axis.text.x  = element_text(angle = 90), axis.title = element_text(size = 28),
          legend.text = element_text(size = 28),
          legend.title = element_text(size = 28)) +
    labs(x='Interaction motifs (top 200 by degree)', y='# of edges (mean degree) - percent', fill = 'Motif source')+
    scale_fill_manual(values = c(rppi_mcolor, rppi_pcolor)) +
    geom_text(mapping = aes(x=motif, y=csum+210, label = sprintf('%s - %s', csum, percent(csum_percent)), angle = 90)) + 
    geom_text(mapping = aes(x=20, y = 3300, label = csum_label), size = 10) +
    geom_text(mapping = aes(x=20, y = 3100, label = csum_label_node), size = 10) +
  outpdf(infile, 'downsampled_ppi_internet_degree_cumsum', width = 28, height = 15)
  ### end cum degree distribution
  # stop()
  print(topndf)
  ###
  # network bit
  ###
  name_parts = paste(strsplit(infile, '_')[[1]][4:7], collapse = '_')
  outname = sprintf('abdb_figures_2019/%s_downsampled_ppi_internet.pdf', name_parts)
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
  node_outfile = 'abdb_outfiles_2019/downsampled_ppi_internet_nodes.csv'
  write.csv(x=nodes, file = node_outfile, row.names = FALSE)
  edge_outfile = 'abdb_outfiles_2019/downsampled_ppi_internet_edges.csv'
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
  # stop()
  col_list = list("Paratope" = '#3399CC', "Epitope"= '#C0C0C0')
  node_attribute[node_attribute=='Paratope'] = rppi_mcolor
  node_attribute[node_attribute=='Epitope'] = rppi_pcolor
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
  rinfile = 'abdb_outfiles_2019/ppi_internet_nodes.csv'
  peinfile = 'abdb_outfiles_2019/downsampled_ppi_internet_nodes.csv'
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
    geom_text(mapping = aes(x=0, y=100, label = scor_label), size=7, hjust=0) + 
    geom_density_2d(mapping = aes(x=score.x, y=score.y), color = my_spectral[3]) +
    labs(x=sprintf('Randomly sampled (degree >= %s)', cutoff), y='Observed node degree') + 
    theme(axis.title = element_text(size=20), panel.background = element_blank())
  outpdf('ppi_ppidownsampled_degree', 'corr')
}


cross_reactivity_density = function() {
  infile = 'abdb_outfiles_2019/sub0.01_ppi_partner_overlap.csv'
  df = read_csv(infile)
  print(df)
  pdf = subset(df, motif_source=='ppimotif')
  rpdf = subset(df, motif_source =='dsppimotif')
  edf = subset(df, motif_source=='ppimotifpartner')
  redf = subset(df, motif_source =='dsppimotifpartner')
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
  ksdf = tibble('motif_source' = c('ppimotif', 'ppimotifpartner'), 'kstest' = c(pks_label, eks_label))
  print(ksdf)
  df[df=='ppimotif'] = 'PPI motif'
  df[df=='ppimotifpartner'] = 'PPI motif partner'
  df[df=='dsppimotif'] = 'Downsampled PPI motif'
  df[df=='dsppimotifpartner'] = 'Downsampled PPI motif partner'
  df$motif_source  = factor(df$motif_source, levels = c('PPI motif', 'PPI motif partner', 'Downsampled PPI motif', 'Downsampled PPI motif partner'))
  meds  = df %>% group_by(motif_source) %>% summarise(med = median(percent_overlap))
  print(df)
  med_labels = sprintf('Median: %s', percent(meds$med/100))
  print(meds)
  #Function to transform data to y positions
  break_point = 0.002
  trans <- function(x){pmin(x,break_point) + 0.0005*pmax(x-break_point,0)}
  yticks= c(0.001, 0.002, 0.7, 2)
  yticks_labels = c('0.001','\\ \\', '0.7', '2')
  print(trans(yticks))
  # stop()
  ggplot(data= df) + 
    geom_density(mapping = aes(x=percent_overlap, y=trans(..density..), fill=motif_source, color = motif_source), color = NA, alpha = 0.6) + 
    geom_text(data = meds, mapping = aes(x=c(2, 2, 14, 14), y = trans(c(2, 1.7, 2, 1.7)), label= med_labels, color=motif_source), hjust=0, show.legend = FALSE) + 
    geom_text(data = ksdf, mapping = aes(x=c(26, 26), y = trans(c(2, 1.7)), label= c(pks_label, eks_label)), hjust =0) + 
    scale_fill_manual(values = c(ppi_mcolor, ppi_pcolor, rppi_mcolor, rppi_pcolor)) +  
    scale_color_manual(values = c(ppi_mcolor, ppi_pcolor, rppi_mcolor, rppi_pcolor)) +  
    scale_y_continuous(limits=c(0,NA), breaks=trans(yticks), labels=yticks_labels) +
    geom_rect(aes(xmin=0, xmax=0.76, ymin=trans(0.05), ymax=trans(0.06)), fill="white") + 
    labs(x= 'Interaction partner overlap (%)', y='Density', fill = 'Motif source', color = 'Motif source')
  outpng('ppi_motif_partner', 'density', width = 20)

}


cross_reactivity_degree_cor = function(){
  rinfile = 'abdb_outfiles_2019/downsampled_ppi_internet_nodes.csv'
  peinfile = 'abdb_outfiles_2019/ppi_internet_nodes.csv'
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
  infile = 'abdb_outfiles_2019/sub0.01_ppi_partner_overlap.csv'
  df = read_csv(infile) %>%.[.$percent_overlap < 100,] # excludes self pairing
  meancross = df %>% group_by(.,motif1, motif_source) %>% summarise('ave' = mean(percent_overlap)) %>% .[order(-.$ave),]
  colnames(meancross)[1] = 'id'
  merge_degree_overlapdf = merge(mdegre_df,meancross, by = 'id')
  print(meancross)
  # stop()
  print(head(merge_degree_overlapdf))
  lvls  = c('ppimotif', 'ppimotifpartner', 'dsppimotif', 'dsppimotifpartner')
  facet_labels = c('PPI motif', 'PPI motif partner', 'Downsampled PPI motif', 'downsampled PPI motif partner')
  names(facet_labels) = lvls
  merge_degree_overlapdf$motif_source = factor(merge_degree_overlapdf$motif_source, levels = lvls)
  for (msource in unique(merge_degree_overlapdf$motif_source)) {merge_degree_overlapdf[]}
  ggplot(data =merge_degree_overlapdf) + 
    geom_point(mapping = aes(x=score, y=ave, color=motif_source), size = 5, alpha=0.5) + 
    facet_wrap(~ motif_source, ncol = 4, labeller = labeller(motif_source=facet_labels), scales = 'free') +
    scale_color_manual(values = c(ppi_mcolor, ppi_pcolor, rppi_mcolor, rppi_pcolor)) + 
    labs(x='Node degree', y = 'Mean interaction partner overlap') + 
    theme(legend.position = 0)
  outpdf('degree', 'ppi_mean_overlap_scatter', width = 12, height = 4)
  
}



# run stuff
ppi_motif_interaction()
# downsampled_ppi_motif_internet()
# degree_correlation()
# cross_reactivity_density()
# cross_reactivity_degree_cor()