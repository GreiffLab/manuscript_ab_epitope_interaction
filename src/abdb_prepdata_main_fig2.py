# import stuff


import pandas as pd
import numpy as np
from mpi4py import MPI
import sys
import numpy as np
import random
import find_files as fifi

#sets df to display all columns
pd.set_option('display.max_column', None)


def victor_cumulative_curve():
    '''
    prepdata for victor's cumulative curve
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged.csv'
    df = pd.read_csv(infile)
    chains = df.abchain.unique()
    data = []
    n = 1000
    pdbids = df.pdbid.unique()
    print(len(pdbids))
    print(df.head())
    #skips some structure so we can see the plot
    skip = 30
    indices = range(0,len(pdbids), skip)
    abag = ['ab_motif', 'ag_motif']
    for motiftype in abag:
        source = motiftype.split('_')[0]
        for run in range(n):
            random.shuffle(pdbids)
            current_motif = []
            for i,pdbid in enumerate(pdbids):
                if i in indices:
                    pdbdf = df[df.pdbid == pdbid]
                    strpdbidx = 's%s' % str(i+2)
                    pdbidx = i+2
                    pdbmotif = getattr(pdbdf,motiftype).tolist()
                    if i == 0:
                        pdbid2 = pdbids[1]
                        pdbdf2 = df[df.pdbid == pdbid2]
                        pdbmotif2 = getattr(pdbdf2,motiftype).tolist()
                        intersection = set(pdbmotif) & set(pdbmotif2)
                        current_motif += list(set(pdbmotif + pdbmotif2))
                        fraction_overlap = float(len(intersection))/len(pdbmotif)
                        num_unique_motif = len(set(current_motif))
                        datum = [pdbidx, strpdbidx,pdbid, fraction_overlap, source, num_unique_motif]
                        print(datum)
                        data.append(datum)
                    elif i == 1:
                        data.append(datum)
                    elif i >=2:
                        intersection = set(pdbmotif) & set(current_motif)
                        fraction_overlap = float(len(intersection))/len(pdbmotif)
                        # print(pdbmotif)
                        # print(intersection)
                        # print(fraction_overlap)
                        current_motif += list(set(pdbmotif))
                        num_unique_motif = len(set(current_motif))
                        datum = [pdbidx, strpdbidx,pdbid, fraction_overlap, source, num_unique_motif]
                        data.append(datum)
                    print('run %s' % run)
    colnames = ['pdbidx', 'strpdbidx', 'pdbname', 'coverage', 'source', 'unique_motif']
    covdf = pd.DataFrame(data, columns=colnames)
    print(covdf)
    outname = 'abdb_outfiles_2019/motif_coverage_skip%s.csv' % skip
    covdf.to_csv(outname, index=False)

def add_motif_len(infile):
    '''
    add motif len to the dataframe
    note: for epitope, gapset = egapset, plen = elen
    :return: 
    '''
    # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_' \
    #          'segments_abshift_abshiftl_epitope_segment_notationx.csv'
    df = pd.read_csv(infile)
    print(df.head())
    df= df.iloc[:]
    lens = []
    for i,row in df.iterrows():
        if 'epitope' in infile:
            gapsets = row.egapset
            length = row.elen
        else:
            gapsets= row.gapset
            length = row.plen
        motif = row.gap_patterns
        if str(gapsets) == 'nan':
            motif_len = 1
            lens.append(motif_len)
        else:
            gapsets = gapsets.split('-')
            gapsets = [item for item in gapsets if int(item) > 0]
            motif_len  = len(gapsets) + length
            lens.append(motif_len)
        print(motif, motif_len)
    df['motif_len'] = lens
    outname = infile.split('.')[0] + '_len.csv'
    df.to_csv(outname, index = False)

def merge_para_epi_notation_files():
    '''
    merge the two files:
    respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_len.csv
    respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_len.csv

    :return:
    '''
    paradf = pd.read_csv('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_' \
                         'segments_abshift_abshiftl_paratope_segment_notationx_len.csv')
    epidf = pd.read_csv('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_'\
                        'segments_abshift_abshiftl_epitope_segment_notationx_len.csv')
    print(epidf.head())
    print(paradf.head())
    epitopes = []
    elens = []
    motif_lens_epi = []
    ag_motifs = []
    agresnumisets = []
    agchains = []
    for i,row in epidf.iterrows():
        epitope = row.epitope
        elen = row.elen
        motif_len_epi = row.motif_len
        epitopes.append(epitope)
        elens.append(elen)
        motif_lens_epi.append(motif_len_epi)
        ag_motif = row.gap_patterns
        ag_motifs.append(ag_motif)
        agresnumiset = row.agresnumiset
        agresnumisets.append(agresnumiset)
        agchains.append(row.agchain)
    paradf['ag_motiflen'] = motif_lens_epi
    paradf['epitope'] = epitopes
    paradf['epitope_len'] = elens
    paradf['ag_motif'] = ag_motifs
    paradf['agresnumiset'] = agresnumisets
    paradf['agchain'] = agchains
    paradf = paradf.rename(columns = {'motif_len': 'ab_motiflen', 'gap_patterns':'ab_motif'})
    outname = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged.csv'
    paradf.to_csv(outname, index=False)

def add_pdb_resolution():
    '''
    add resolution for each pdb structure
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged.csv'
    df = pd.read_csv(infile)
    print(df.head())
    resolutions = []
    data = []
    pdbids = df.pdbid.unique()
    for pdbid in pdbids:
        pdbfile = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin/%s.pdb' % pdbid
        contents = open(pdbfile).readlines()
        resolution = contents[2].split()[-1].strip()
        print(resolution)
        data.append([pdbid, resolution])
    outdf = pd.DataFrame(data, columns= ['pdbid', 'resolution'])
    print(outdf)
    outname = 'abdb_outfiles_2019/abdb_resolution.csv'
    outdf.to_csv(outname, index=False)




# run stuff
# victor_cumulative_curve()
add_motif_len('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx.csv')
add_motif_len('abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx.csv')
merge_para_epi_notation_files()
add_pdb_resolution()




