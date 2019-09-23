# notationx merging and max len are done in fram: /cluster/projects/nn9603k/icml_seq/mat/abdb/3did/take2/merged_files
# import stuff
import pandas as pd
import sys

# set display to max
pd.set_option('display.max_column', None)


# one to three letter amino acid dict
aafile = '../datasets/amino_acids/the_twenty.txt'
aacontent = open(aafile).read().splitlines()
aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)


def crosscheck_abdb_3did():
    '''
    check pdb id overlap between abdb and 3did
    output residue pairs from 3did
    :return:
    '''
    abdbfile = 'abdb_outfiles_2019/respairs_absort_cutoff5.csv'
    abdbpdbid = pd.read_csv(abdbfile).pdbid.unique().tolist()
    abdbpdbid = [item.split('_')[0] for item in abdbpdbid]
    tdidfile = '../datasets/3did/3did_flat.txt'
    lines = open(tdidfile).read().splitlines()
    print(len(lines))
    newcontents = []
    tdparts = lines[1].split('\t')
    currentpdbid = tdparts[1].upper()
    current_chain1 = tdparts[2][0]
    current_chain2 = tdparts[3][0]
    idparts = lines[0].split('\t')
    current_domain1 = idparts[1]
    current_domain2 = idparts[2]
    lines = [item for item in lines if '//' not in item]
    alphabets = list(string.ascii_uppercase)
    alphabets = dict([item,i+1] for i,item in enumerate(alphabets))
    for line in lines[2:]:
        if '#=3D' not in line and '#=ID'not in line:
            parts = [currentpdbid, current_domain1, current_domain2, current_chain1, current_chain2] + \
                    [item.strip() for item in line.split('\t')[:4]]
            #handle insertion
            resnums = parts[-2:]
            parts_inserts = []
            for resnum in resnums:
                try:
                    resnum = int(resnum)
                    parts_inserts += [resnum, 0]
                except:
                    insertchar = resnum[-1].upper()
                    insertion = alphabets[insertchar]
                    resnum = int(resnum[:-1])
                    parts_inserts += [resnum, insertion]
            if current_chain1 == current_chain2:
                parts = parts[:-2] + parts_inserts + ['intradomain']
            else:
                parts = parts[:-2] + parts_inserts +  ['interdomain']
        elif '#=3D' in line:
            tdparts = line.split('\t')
            currentpdbid = tdparts[1].upper()
            current_chain1 = tdparts[2][0]
            current_chain2 = tdparts[3][0]
        elif '#=ID' in line:
            idparts = line.split('\t')
            current_domain1 = idparts[1]
            current_domain2 = idparts[2]
        if parts[0] not in abdbpdbid:
            # filter out abdb data
            pdbchainpair = parts[0] + '_' + parts[3] + parts[4]
            resnumi1  = str(parts[7]) + '_' + str(parts[8])
            resnumi2  = str(parts[9]) + '_' + str(parts[10])
            parts = parts + [pdbchainpair, resnumi1, resnumi2]
            print(len(parts))
            print(parts)
            # sys.exit()
            newcontents.append(parts)
        # newcontents.append(parts)
    colnames = ['pdbid', 'domain1', 'domain2', 'chain1', 'chain2', 'res1', 'res2', \
                'resnum1', 'ins1', 'resnum2', 'ins2', 'intertype', 'pdbchainpair', 'resnumi1', 'resnumi2']
    outdf = pd.DataFrame(newcontents, columns=colnames)
    print(outdf.info())
    tdidpdbid = outdf.pdbid.unique().tolist()
    overlap = list(set(tdidpdbid) & set(abdbpdbid))
    print('overlap with abdb: %s' % len(overlap))
    outname = 'abdb_outfiles_2019/threedid_no_ig.csv'
    outdf.to_csv(outname, index=False)
    return outname




def filter_ppi_data():
    '''
    filters ppi data (3did) for max_gap < 8 and motif len <= 300
    this ensure that the ppi data does not drift too much from abdbd data
    :return:
    '''
    infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged.csv'
    df = pd.read_csv(infile)
    print(df.shape)
    max_gap = 7
    max_len = 300
    df = df[(df.max_gap <= max_gap) & (df.motif_len <= max_len)]
    outname = infile.split('.')[0] + '_maxgap%s_maxlen%s.csv' % (max_gap, max_len)
    df.to_csv(outname, index=False)

def get_ppi_residues():
    '''
    get ppi three letter residue names
    :return:
    '''
    infile  = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300.csv'
    df = pd.read_csv(infile).iloc[:]
    print(df.shape)
    df = df[(df.max_gap < 8)]
    print(df.shape)
    data = []
    aadictr = dict([(aadict[item],item) for item in aadict.keys()])
    for i, row in df.iterrows():
        sequence = row.sequence
        residues = list(sequence)
        residuestl = [aadictr[item] for item in residues]
        for residue in residuestl:
            datum = [row.pdbchainpair,residue, row.intertype]
            data.append(datum)
    colnames = ['pdbchainpair', 'residue', 'domain_type']
    outdf = pd.DataFrame(data, columns=colnames)
    outname = infile.split('.')[0] + '_three_letter_residue.csv'
    outdf.to_csv(outname, index=False)
    print(outdf.head())


def find_immunofam():
    '''
    find immunoglobulin and Ig-like domain from pfam
    :return:
    '''
    infile = '../datasets/pfam/pdb_pfam_mapping.txt'
    df = pd.read_csv(infile, delimiter='\t')
    keywords  = ['Immunoglobulin', 'Ig-like']
    newcontent = []
    for i, row in df.iterrows():
        for keyword in keywords:
            if keyword in row.PFAM_desc:
                newcontent.append(row)
    igdf = pd.DataFrame(newcontent, columns=df.columns)
    igfam = igdf.PFAM_Name.unique()
    # print(igfam, len(igfam))
    return igfam

def filter_iglike():
    '''
    filter out immunoglobulin and iglike structures
    output no_iglike file
    :return:
    '''
    igfam = find_immunofam()
    print(igfam)
    infile = 'abdb_outfiles_2019/threedid_no_ig.csv'
    df = pd.read_csv(infile)
    print('with iglike: %s' % df.shape[0])
    df2 = df[(~df.domain1.isin(igfam)) & (~df.domain2.isin(igfam))]
    print('no iglike: %s' % df2.shape[0])
    outname = infile .split('.')[0] + 'like.csv'
    df2.to_csv(outname, index=False)

def get_paired_motifs():
    '''
    get paired motif from noiglike and res/size filtered ppi motifs
    :return:
    '''
    infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300.csv'
    df = pd.read_csv(infile)
    print(df.head())
    print(df.shape)
    print(len(df.gap_pattern.unique()))
    ddf = df[df.duplicated(subset = 'pdbchainpair')]
    print(ddf.shape)
    print(df.head())
    data = []
    for pdbchainpair in ddf.pdbchainpair.tolist()[:]:
        cpdf = df[df.pdbchainpair == pdbchainpair]
        datum = cpdf.iloc[0].tolist() + cpdf.iloc[1].tolist()
        data.append(datum)
    colnames = [item+'1' for item in cpdf.columns] + [item+'2' for item in cpdf.columns]
    outdf = pd.DataFrame(data, columns=colnames)
    print(outdf.head())
    outname = infile.split('.')[0] + '_paired.csv'
    outdf.to_csv(outname, index=False)



# run stuff
# crosscheck_abdb_3did()
# find_immunofam()
# filter_iglike()
# filter_ppi_data()
# get_ppi_residues()
get_paired_motifs()
