# import stuff
import jellyfish
import pandas as pd
import sys

def get_levenshtein_segments():
    '''
    get levenshtein distance per segment
    :return:
    '''
    infile = 'abdb_outfiles_2019/abdb_segment_absequence_full_vgene_imgt_vgene.csv'
    df = pd.read_csv(infile)
    print(df.info())
    data = []
    for segment in df.segment.unique():
        segdf = df[df.segment == segment]
        print(segment)
        print(segdf.shape)
        counter = 0
        for i, row in segdf.iterrows():
            counter += 1
            print(counter)
            seq1 = row.segment_seq
            pdbid = row.pdbid
            for i2, row2 in segdf.iterrows():
                pdbid2 = row2.pdbid
                if pdbid != pdbid2:
                    seq2 = row2.segment_seq
                    ld = jellyfish.levenshtein_distance(seq1,seq2)
                    datum = [pdbid, pdbid2, segment, seq1, seq2,ld]
                    data.append(datum)
    colnames = ['pdbid1', 'pdbid2', 'segment', 'seq1', 'seq2', 'ld']
    lddf = pd.DataFrame(data, columns=colnames)
    print(lddf.head())
    outname = infile[:-4] + '_ld.csv'
    print(outname)
    lddf.to_csv(outname, index=False)


def add_full_sequence_ab_ag():
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
    df = pd.read_csv(infile)
    pdbdir = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin/'
    aa3to1 = {
        'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M',
        'ILE': 'I', 'LEU': 'L', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
        'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H',
        'CYS': 'C', 'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G',
        'MSE': 'M',
    }
    lseqs = []
    hseqs = []
    aseqs = []
    for i,row in df.iterrows():
        pdbid = row.pdbid
        pdb_file = pdbdir + pdbid + '.pdb'
        print(pdb_file)
        ca_pattern = re.compile(
            "^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
        filename = os.path.basename(pdb_file).split('.')[0]
        chain_dict = dict()
        chain_list = []

        fp = open(pdb_file, 'rU')
        for line in fp.read().splitlines():
            if line.startswith("ENDMDL"):
                break
            match_list = ca_pattern.findall(line)
            if match_list:
                resn = match_list[0][0] + match_list[0][2]
                chain = match_list[0][1] + match_list[0][3]
                if chain in chain_dict:
                    chain_dict[chain] += aa3to1[resn]
                else:
                    chain_dict[chain] = aa3to1[resn]
                    chain_list.append(chain)
        ab_chain = ['L', 'H']
        ag_chain = [item for item in chain_list if item not in ab_chain]
        lseqs.append(chain_dict['L'])
        hseqs.append(chain_dict['H'])
        aseqs.append(chain_dict[ag_chain[0]])
    outdf = pd.DataFrame()
    outdf['pdbid'] = df.pdbid
    outdf['region'] = df.segment
    outdf['l_sequence'] = lseqs
    outdf['h_sequence'] = hseqs
    outdf['paratope'] = df.paratope
    outdf['a_sequence'] = aseqs
    outdf['epitope']  = df.epitope
    print(outdf.head())
    outname = 'abdb_outfiles_2019/' + 'heavy_light_ag_aaseq.csv'
    outdf.to_csv(outname , index=False)


def get_levenshtein_agseq():
    '''
    get levenshtein distance per antigen
    :return:
    '''
    infile = 'abdb_outfiles_2019/heavy_light_ag_aaseq.csv'
    df = pd.read_csv(infile).iloc[:]
    print(df.info())
    data = []
    for i, pdbid in enumerate(df.pdbid.unique()):
        pdbdf = df[df.pdbid == pdbid]
        agseq1 = pdbdf.iloc[0].a_sequence
        print('computing %s #%s' % (pdbid, i) )
        for pdbid2 in df.pdbid.unique():
            if pdbid2 != pdbid:
                pdbdf2 = df[df.pdbid == pdbid2]
                agseq2 = pdbdf2.iloc[0].a_sequence
                ld = jellyfish.levenshtein_distance(agseq1,agseq2)
                # print(ld)
                datum = [pdbid, pdbid2, agseq1, agseq2, ld]
                data.append(datum)
    colnames = ['pdbid1', 'pdbid2', 'agseq1', 'agseq2', 'ld']
    lddf = pd.DataFrame(data, columns=colnames)
    outname = infile[:-4] + '_antigen_full_ld.csv'
    print(outname)
    lddf.to_csv(outname, index=False)


def get_levenshtein_epitopeseq():
    '''
    get levenshtein distance per antigen
    :return:
    '''
    infile = 'abdb_outfiles_2019/heavy_light_ag_aaseq.csv'
    df = pd.read_csv(infile).iloc[:]
    print(df.info())
    data = []
    for i, row in df.iterrows():
        pdbid = row.pdbid
        epitopeseq1 = row.epitope
        for i2, row2 in df.iterrows():
            pdbid2 = row2.pdbid
            if pdbid2 != pdbid:
                epitopeseq2 = row2.epitope
                ld = jellyfish.levenshtein_distance(epitopeseq1,epitopeseq2)
                datum = [pdbid, pdbid2, epitopeseq1, epitopeseq2, ld]
                data.append(datum)
    colnames = ['pdbid1', 'pdbid2', 'epitopeseq1', 'epitopeseq2', 'ld']
    lddf = pd.DataFrame(data, columns=colnames)
    print(lddf.head())
    outname = infile[:-4] + '_antigen_epitope_ld.csv'
    print(outname)
    lddf.to_csv(outname, index=False)


def get_levenshtein_segments_epitope():
    '''
    get levenshtein distance per segment
    :return:
    '''
    infile = 'abdb_outfiles_2019/heavy_light_ag_aaseq.csv'
    df = pd.read_csv(infile).iloc[:]
    print(df.info())
    df = df.dropna(subset = ['epitope'])
    data = []
    for segment in df.region.unique():
        segdf = df[df.region == segment]
        print(segment)
        print(segdf.shape)
        counter = 0
        for i, row in segdf.iterrows():
            counter += 1
            # print(counter)
            print('seq1 %s' % row.epitope)
            seq1 = row.epitope
            pdbid = row.pdbid
            for i2, row2 in segdf.iterrows():
                pdbid2 = row2.pdbid
                if pdbid != pdbid2:
                    print('seq2 %s' % row2.epitope)
                    seq2 = row2.epitope
                    ld = jellyfish.levenshtein_distance(seq1,seq2)
                    datum = [pdbid, pdbid2, segment, seq1, seq2,ld]
                    data.append(datum)
    colnames = ['pdbid1', 'pdbid2', 'region', 'epitope1', 'epitope2', 'ld']
    lddf = pd.DataFrame(data, columns=colnames)
    print(lddf.head())
    outname = infile[:-4] + '_antigen_epitope_ld.csv'
    print(outname)
    lddf.to_csv(outname, index=False)


# run stuff
get_levenshtein_segments()
# add_full_sequence_ab_ag()
# get_levenshtein_agseq()
# get_levenshtein_epitopeseq()
# get_levenshtein_segments_epitope()