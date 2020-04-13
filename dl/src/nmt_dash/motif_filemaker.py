# make motif input file with the following format
# motif \t motif
import pandas as pd
import sys

infile = '/Users/rahmadakbar/greifflab/aims/aimugen/src/abdb_outfiles/respairs_segment_notationx_merged_angle_len.csv'
df = pd.read_csv(infile)
# print(df.head)
df = df.dropna(subset=['paratope', 'epitope'])
outcontent = ''
for i, row in df.iterrows():
    agmotif = row.ag_gap_pattern
    abmotif = row.ab_gap_pattern
    # print(agmotif, abmotif)
    abparts = abmotif.split('X')
    agparts  = agmotif.split('X')
    abcontents = []
    for abpart in abparts[:-1]:
        if abpart == '':
            abcontents.append('X')
        else:
            abcontents.append(abpart)
            abcontents.append('X')
    agcontents = []
    for agpart in agparts[:-1]:
        if agpart == '':
            agcontents.append('X')
        else:
            agcontents.append(agpart)
            agcontents.append('X')

    # print(abcontents)
    # print(agcontents)
    content = '\t'.join([' '.join(agcontents), ' '.join(abcontents)]) + '\n'
    # print([content])
    outcontent += content
outpath = '../../dataset/motif_epipara.tsv'
outfile = open(outpath, 'w')
outfile.write(outcontent)
