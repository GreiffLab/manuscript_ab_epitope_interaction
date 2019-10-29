load ../datasets/NR_LH_Protein_Martin/2DQC_1.pdb
select me, chain H and resi 95+96
color bluewhite, all
color blue, chain H
color cyan, chain L
color salmon, me
disable me
zoom me, 10
ray 2000
save abdb_figures_2019/2DQC_1_0.png
bg_color white