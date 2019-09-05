# import stuff
import pandas as pd
import sys

pd.set_option('display.max_column', None)

def get_gap_residues_epitope():
    '''
    get residues comprising gaps
    :return:
    '''
    pdbdir = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin/'
    infile = \
        'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_len' \
        '.csv'
    df = pd.read_csv(infile).iloc[:]
    df = df[df.gapstrstatus != 'continuous']
    print(df.head())
    data = []
    for i, row in df.iterrows():
        pdbid = row.pdbid
        motif = row.gap_patterns
        chain = row.agchain
        gaps = row.egapset.split('-')
        agresnums = row.agresnumiset.split('-')
        gap_agresnums = [item for item in zip(gaps, agresnums[:-1], agresnums[1:]) if int(item[0]) != 0]
        gapresnums  = []
        for gap_agresnum in gap_agresnums:
            gap = int(gap_agresnum[0])
            start = min([int(item) for item in gap_agresnum[1:]])
            gapresnum = [start + i for i in range(1, gap +1)]
            print(agresnums)
            print(gapresnum)
            gapresnums += gapresnum
        pdbfile = pdbdir + pdbid + '.pdb'
        lines = open(pdbfile).readlines()
        for line in lines:
            if line.startswith('ATOM'):
                line_chain = line[21]
                if line_chain == chain:
                    resnum = int(line[22:26].strip())
                    if resnum in gapresnums:
                        residue = line[17:20].strip()
                        print(residue)
                        datum = [residue, row.segment, row.abchain]
                        print(datum)
                        data.append(datum)
    colnames = ['residue', 'segment', 'abchain']
    outdf = pd.DataFrame(data, columns=colnames)
    outname = infile.split('.')[0] + '_gap_residue.csv'
    outdf.to_csv(outname, index=False)



def get_gap_residues_paratope():
    '''
    get residues comprising gaps
    :return:
    '''
    pdbdir = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin/'
    infile = \
        'abdb_outfiles_2019' \
        '/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx_len' \
        '.csv'
    df = pd.read_csv(infile).iloc[:]
    df = df[df.gapstrstatus != 'continuous']
    print(df.head())
    data = []
    for i, row in df.iterrows():
        pdbid = row.pdbid
        motif = row.gap_patterns
        chain = row.abchain
        gaps = row.gapset.split('-')
        abresnums = row.abresnumiset.split('-')
        gap_abresnums = [item for item in zip(gaps, abresnums[:-1], abresnums[1:]) if int(item[0]) != 0]
        gapresnums  = []
        for gap_abresnum in gap_abresnums:
            gap = int(gap_abresnum[0])
            print(abresnums)
            try:
                start = min([int(item) for item in gap_abresnum[1:]])
            except:
                start = int(gap_abresnum[1][:-1])
            gapresnum = [start + i for i in range(1, gap +1)]
            print(gapresnum)
            gapresnums += gapresnum
        pdbfile = pdbdir + pdbid + '.pdb'
        lines = open(pdbfile).readlines()
        for line in lines:
            if line.startswith('ATOM'):
                line_chain = line[21]
                if line_chain == chain:
                    resnum = int(line[22:26].strip())
                    if resnum in gapresnums:
                        residue = line[17:20].strip()
                        print(residue)
                        datum = [residue, row.segment, row.abchain]
                        print(datum)
                        data.append(datum)
    colnames = ['residue', 'segment', 'abchain']
    outdf = pd.DataFrame(data, columns=colnames)
    outname = infile.split('.')[0] + '_gap_residue.csv'
    outdf.to_csv(outname, index=False)



# run stuff
# get_gap_residues()
get_gap_residues_paratope()