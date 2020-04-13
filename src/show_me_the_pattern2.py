#takes position as the first argument and segment as the second argument

# import stuff

import os
import sys
import pandas as pd
import random

#pattern = '0-0'
pattern = sys.argv[1]
segment = sys.argv[2]
df = pd.read_csv('abdb_outfiles/respairs_paratope_segment.csv')
print(df.head())
pdf = df[(df.gapset == pattern) & (df.segment == segment)]
print(pdf)
random_pdbid = random.choice(pdf.pdbid.values)
pdbdf = pdf[pdf.pdbid == random_pdbid]
resnums = '+'.join(pdbdf.abresnumiset.values[0].split('-'))
chain = pdbdf.abchain.values[0]
random_pdbid_path = '../datasets/NR_LH_Protein_Martin/%s.pdb' % random_pdbid
pmlcontents = []
#pmlcontent = 'load %s\nselect me, chain H and resi %s\nshow spheres, me\nutil.cbc\nbg_color white\ndisable me'%(
#    random_pdbid_path,resnums)
pmlcontents.append('load %s' % random_pdbid_path)
pmlcontents.append('select me, chain %s and resi %s' % (chain,resnums))
pmlcontents.append('color bluewhite, all')
pmlcontents.append('color blue, chain H')
pmlcontents.append('color cyan, chain L')
#pmlcontents.append('show dots, me')
pmlcontents.append('color salmon, me')
pmlcontents.append('disable me')
pmlcontents.append('zoom me, 10')
pmlcontents.append('ray 2000')
outimg = 'abdb_figures_2019/%s_%s.png' % (random_pdbid, sys.argv[1])
pmlcontents.append('save %s' % outimg)
pmlcontents.append('bg_color white')
#pmlcontents.append('quit')
pmlcontent = '\n'.join(pmlcontents)
print(pmlcontent)
outfile = open('showme.pml', 'w')
outfile.write(pmlcontent)
outfile.close()
os.system('pymol showme.pml')
os.system('open %s' % outimg)
