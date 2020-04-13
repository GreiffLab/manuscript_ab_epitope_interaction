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

