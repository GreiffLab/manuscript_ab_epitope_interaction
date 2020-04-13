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

### merged from prepdata_main_figure4.py


# import stuff
import sys
import pandas as pd


pd.set_option('display.max_column', None)

def get_context_pattern_top3(infile):
    '''
    gets context from the union of top 3 epitope and paratope motifs.
    motifs: XXX, X1X, XX, X2X
    no longer top 3, we do 4 now
    :return:
    '''
    df = pd.read_csv(infile)
    print(df.shape)
    if 'epitope' in infile:
        df['paratope'] = df.epitope
    dfx = df.dropna(subset=['gap_patterns', 'paratope'])
    print(dfx.shape)
    patterns = ['XXX', 'X1X', 'X2X', 'XX']
    print(patterns)
    ptag = '_'.join(patterns)
    dfxs = []
    for pattern in patterns:
        dfp = dfx[df.gap_patterns == pattern]
        dfxs.append(dfp)
    dfxc = pd.concat(dfxs)
    print(dfxc.shape)
    edges = []
    edges_next = []
    edges_next_single = []
    edges_prev_next_single = []
    for i,row in dfxc.iterrows():
        paratope = row.paratope
        # print(row)
        motif = row.gap_patterns
        lenp = len(paratope)
        for i2, res in enumerate(paratope):
            next_i = -(lenp-i2-1)
            next_res = paratope[next_i:]
            prev_res = paratope[:i2]
            if i2 == 0:
                prev_res = '-'
            if next_res == paratope: # check for last residue
                context = prev_res
                context_next_single = '-'
                context_prev_next_single = prev_res[-1] + '-'
            else:
                context = prev_res + next_res
                context_next = next_res
                context_next_single = next_res[0]
                context_prev_next_single = prev_res[-1] + next_res[0]
                if res == context_next_single:
                    context_next_single += "'"
            edge = (res, context, motif)
            edge_next = (res, context_next, motif)
            edge_next_single = (res, context_next_single, motif)
            edge_prev_next_single = (res, context_prev_next_single, motif)
            edges.append(edge)
            edges_next.append(edge_next)
            if i2 != lenp-1:
                edges_next_single.append(edge_next_single)
            edges_prev_next_single.append(edge_prev_next_single)
            print(paratope, res, context_prev_next_single, context_next_single)
        print('%s/%s: %s' % (i, dfxc.shape[0],infile))
    print(edges_next_single)
    colnames = ['source', 'target', 'motif']
    edge_nextdf = pd.DataFrame(edges_next_single, columns=colnames)
    print(edge_nextdf.head())
    pairdict1 = {}
    sources = []
    targets = []
    for i, row in edge_nextdf.iterrows():
        key1 = row.source + row.target
        key2 = row.target + row.source
        val = row.source + row.target
        if key1 not in pairdict1:
            pairdict1[key1] = val
        if key2 not in pairdict1:
            pairdict1[key2] = val
        if key2 in pairdict1:
            source = pairdict1[key2][0]
            target = pairdict1[key2][1]
            sources.append(source)
            targets.append(target)
        else:
            sources.append(row.source)
            targets.append(row.target)
    print(len(sources))
    print(edge_nextdf.shape)
    edge_nextdf['source2'] = sources
    edge_nextdf['target2'] = targets
    print(edge_nextdf)
    outname = infile.split('.')[0] + '_%s_edge_next.csv' % ptag
    edge_nextdf.to_csv(outname, index=False)



## run stuff
# epiinfile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments' \
#              '_abshift_abshiftl_epitope_segment_notationx.csv'

parainfile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments' \
             '_abshift_abshiftl_paratope_segment_notationx.csv'
# get_context_pattern_top3(epiinfile)
get_context_pattern_top3(parainfile)





