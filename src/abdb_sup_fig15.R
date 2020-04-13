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

theme_set(themeakbar())
my_spectral <- colorRampPalette(brewer.pal(8,'Spectral'))(14)
wgrandbudapest = wes_palette('GrandBudapest2')
print(wgrandbudapest)
my_grandbudapest <- colorRampPalette(wgrandbudapest)(14)
print(my_grandbudapest)
abcolor = '#6464FF'
abcolor2 = '#1010FF'
abcolor3 = '#3399CC'
agcolor = '#C0C0C0'
blues = colorRampPalette(brewer.pal(8,'Blues'))(14)[2:10]
oranges = brewer.pal(6,'Oranges')
greys = brewer.pal(9,'Greys')
print(greys)
agcolor2 = greys[6]
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

residue_contact = function(infile, ncol){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_residue_contact_odds.csv'
  xlabel = 'Paratope'
  ylabel = 'Epitope'
  height = 14
  width = 25
  axistitlesize = 45 
  legendtitlesize = 45
  legendtextsize = 45
  striptextsize = 40
  if (grepl('threedid',infile)){ 
    xlabel = 'Binding partner' 
    ylabel = 'Residue'
    height = 7
    width = 12.5
    axistitlesize = 25
    legendtitlesize = 25
    legendtextsize = 25
    striptextsize = 25
    }
  df = read_csv(infile)
  sdf = spread(df, key = agres, value = logodd)
  residues = factor(unique(df$abres))
  print(residues)
  reslist = list(charged = c('ARG', 'LYS', 'ASP', 'GLU', 'HIS'), polar = c('GLN', 'ASN', 'SER', 'THR'),
                 amphiphatic = c('PHE','TRP', 'TYR'), nonpolar = c('CYS','ALA', 'ILE', 'LEU', 'MET', 'VAL', 'PRO', 'GLY') )
  df$abres = factor(df$abres, levels = c(reslist$charged, reslist$polar, reslist$amphiphatic, reslist$nonpolar))
  df$agres = factor(df$agres, levels = c(reslist$charged, reslist$polar, reslist$amphiphatic, reslist$nonpolar))
  axiscolor= c(rep('red', length(reslist$charged)), rep('salmon',length(reslist$polar)), rep('orange',length(reslist$amphiphatic)),
               rep('blue', length(reslist$nonpolar)))
  print(axiscolor)
  print(length(reslist$charged))
  collabels  = round(c(min(df$logodd, na.rm = TRUE)+1, 0, max(df$logodd, na.rm = TRUE)-1))
  print(collabels)
  ggplot(data = df) + 
    geom_tile(data= subset(df, !is.na(logodd)), mapping = aes(x = abres, y=agres, fill = logodd), stat = 'identity') + 
    geom_text(mapping = aes(x= abres, y=agres, label = round(logodd))) +
    facet_wrap(~ region, ncol = ncol) +
    theme(axis.text.x = element_text(angle = 90, color = axiscolor),
          axis.text.y = element_text(color=axiscolor),
          panel.background = element_blank(),
          axis.title = element_text(hjust = 0.0, size = axistitlesize),
          legend.text = element_text(size = legendtextsize),
          legend.title = element_text(size = legendtitlesize),
          strip.text = element_text(size = striptextsize)) + 
    scale_fill_gradient2(low = my_spectral[14], mid='white', high = my_spectral[1], breaks = collabels,
                          labels = paste("", collabels),
                          guide = guide_legend(
                            direction = "horizontal",
                            # title.position = "top",
                            # label.position = "bottom",
                            label.hjust = 0.5,
                            label.vjust = 1,)) + 
    labs(x = xlabel , y= ylabel, fill = 'Log2 odds')
  
  outpdf(infile, 'heatmap', height =height, width = width)
}


residue_contact_comp = function(infile, ncol){
  # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_residue_contact_odds.csv'
  xlabel = 'Paratope'
  ylabel = 'Epitope'
  if (grepl('threedid',infile)){ 
    xlabel = 'Binding partner' 
    ylabel = 'Residue'}
  df = read_csv(infile) %>% mutate(respair = paste0(abres,'_',agres))
  print(df)
  compdf  = df %>% group_by(abres, agres) %>% summarise(logoddsum = sum(logodd, na.rm = TRUE))
  print(compdf)
  df = compdf
  residues = factor(unique(df$abres))
  print(residues)
  # reslist = list(charged = c('ARG', 'LYS', 'ASP', 'GLU'), polar = c('GLN', 'ASN', 'HIS', 'SER', 'THR', 'CYS'),
  #                amphiphatic = c('TRP', 'TYR'),hydrophobic = c('ALA', 'ILE', 'LEU', 'MET', 'PHE', 'VAL', 'PRO', 'GLY') )
  # df$abres = factor(df$abres, levels = c(reslist$charged, reslist$polar, reslist$amphiphatic, reslist$hydrophobic))
  # df$agres = factor(df$agres, levels = c(reslist$charged, reslist$polar, reslist$amphiphatic, reslist$hydrophobic))
  # axiscolor= c(rep('red', length(reslist$charged)), rep('salmon',length(reslist$polar)), rep('orange',length(reslist$amphiphatic)),
  #              rep('blue', length(reslist$hydrophobic)))
  reslist = list(charged = c('ARG', 'LYS', 'ASP', 'GLU', 'HIS'), polar = c('GLN', 'ASN', 'SER', 'THR'),
                 amphiphatic = c('PHE','TRP', 'TYR'), nonpolar = c('CYS','ALA', 'ILE', 'LEU', 'MET', 'VAL', 'PRO', 'GLY') )
  df$abres = factor(df$abres, levels = c(reslist$charged, reslist$polar, reslist$amphiphatic, reslist$nonpolar))
  df$agres = factor(df$agres, levels = c(reslist$charged, reslist$polar, reslist$amphiphatic, reslist$nonpolar))
  axiscolor= c(rep('red', length(reslist$charged)), rep('salmon',length(reslist$polar)), rep('orange',length(reslist$amphiphatic)),
               rep('blue', length(reslist$nonpolar)))
  print(axiscolor)
  print(length(reslist$charged))
  collabels  = round(c(min(df$logoddsum)+1, 0, max(df$logoddsum)-1))
  print(collabels)
  # collabels = c(-2,0,1,8)
  # stop()
  ggplot(data = df) + 
    geom_tile(mapping = aes(x = abres, y=agres, fill = logoddsum), stat = 'identity') + 
    geom_text(mapping = aes(x= abres, y=agres, label = round(logoddsum))) +
    # facet_wrap(~ region, ncol = ncol) + 
    theme(axis.text.x = element_text(angle = 90, color = axiscolor),
          axis.text.y = element_text(color=axiscolor),
          panel.background = element_blank(),
          axis.title = element_text(hjust = 0.0, size = 25),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25)) + 
    scale_fill_gradient2(low = my_spectral[14], mid='white', high = my_spectral[1], breaks = collabels,
                          labels = paste("", collabels),
                          guide = guide_legend(
                            direction = "horizontal",
                            # title.position = "top",
                            # label.position = "bottom",
                            label.hjust = 0.5,
                            label.vjust = 1,)) + 
    labs(x = xlabel , y= ylabel , fill = 'Log2 odds')
  
  outpdf(infile, 'heatmap_comp', height =7.3, width = 7)
}



# run stuff
# residue_contact('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_residue_contact_odds.csv', 6)
residue_contact('abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_residue_contact_odds.csv', 2)
# residue_contact_comp('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_residue_contact_odds.csv', 6)
# residue_contact_comp('abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_residue_contact_odds.csv', 2)






