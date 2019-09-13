rsync -av /Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin datasets
rsync -av /Users/rahmadakbar/greifflab/aims/aimugen/datasets/3did/3did_flat.txt.zip datasets/3did
rsync -av /Users/rahmadakbar/greifflab/aims/aimugen/src/abdb_sup_fig*.R src
rsync -av /Users/rahmadakbar/greifflab/aims/aimugen/src/abdb_main_fig*.R src	
rsync -av /Users/rahmadakbar/greifflab/aims/aimugen/src/abdb_prepdata_main_fig*.py src 	
rsync -av /Users/rahmadakbar/greifflab/aims/aimugen/src/abdb_prepdata_sup_fig*.py src 	
rsync -av /Users/rahmadakbar/Google\ Drive/Manuscripts/2018/abdb_manuscript/figures/abdb_main_figure_*.png figures
rsync -av /Users/rahmadakbar/Google\ Drive/Manuscripts/2018/abdb_manuscript/figures/abdb_sup_figure_*.png figures
rsync -av /Users/rahmadakbar/greifflab/aims/aimugen/src/abdb_outfiles_2019/*.csv datasets/preprocessed/
rm src/*conflicted* 	


