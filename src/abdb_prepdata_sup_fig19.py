# import stuf

import sys
import pandas as pd
import os
from abdb_prepdata_main_fig6 import make_tabsep_ppi

def branch_dataset():
    '''
    output sequence dataset for largest paratope branches (XXX)
    :return:
    '''
    nodesfile = 'abdb_outfiles_2019/paratope_epitope_internet_nodes.csv'
    df = pd.read_csv(nodesfile)
    n = 3
    topmotifs = df.id.iloc[:n]
    abagfile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber.csv'
    abagdf = pd.read_csv(abagfile)
    print(abagdf.info())
    for topmotif in topmotifs:
        tdf = abagdf[abagdf.ab_motif == topmotif]
        print(tdf.info())
        datasetdir = '/Users/rahmadakbar/greifflab/aims/aimugen/dl/dataset_%s' % topmotif
        os.system('mkdir %s' % datasetdir)
        outname = '%s/paraepi.tsv' % datasetdir
        outname2 = '%s/epipara.tsv' % datasetdir
        paras = tdf.paratope
        epis = tdf.epitope
        outcontent = '\n'.join(['\t'.join(item) for item in zip(paras, epis)])
        outcontent2 = '\n'.join(['\t'.join(item) for item in zip(epis, paras)])
        outfile  = open(outname, 'w')
        outfile.write(outcontent)
        outfile2 = open(outname2, 'w')
        outfile2.write(outcontent2)



def uniquepair_ppi_dataset():
    '''
    get unique pairs from ppi dataset
    :return:
    '''
    infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv'
    df = pd.read_csv(infile)
    pairs = []
    for i,row in df.iterrows():
        seq1 = row.sequence1
        seq2 = row.sequence2
        pair = '-'.join([seq1,seq2])
        pairs.append(pair)
    print(pairs)
    df['seqpair'] = pairs
    print(df.shape)
    df = df.drop_duplicates(subset='seqpair')
    print(df.shape)
    outfile = infile.split('.')[0] + '_uniqueseqpair.csv'
    df.to_csv(outfile, index=False)
    outseq1 = '../dl/dataset_ppi_uniqueseqpair/paraepi.tsv'
    outseq2 = '../dl/dataset_ppi_uniqueseqpair/epipara.tsv'
    outcontent1 = '\n'.join(['\t'.join([seq1,seq2]) for seq1,seq2 in zip(df.sequence1, df.sequence2)])
    outcontent2 = '\n'.join(['\t'.join([seq2,seq1]) for seq1,seq2 in zip(df.sequence1, df.sequence2)])
    print(outcontent1)
    print(outcontent2)
    outfile1 = open(outseq1, 'w')
    outfile1.write(outcontent1)
    outfile2 = open(outseq2, 'w')
    outfile2.write(outcontent2)

# run stuff
# branch_dataset()
uniquepair_ppi_dataset()