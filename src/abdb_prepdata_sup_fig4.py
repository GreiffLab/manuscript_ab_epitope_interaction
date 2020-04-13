# prepdata for sup fig 4

# import stuff
import sys
import os
from find_files import find_files
import pandas as pd
import re
from Bio import PDB
import jellyfish


#sets df to display all columns
pd.set_option('display.max_column', None)

# create outdir
outpath = 'abdb_outfiles_2019'

# one to three letter amino acid dict
aafile = '../datasets/amino_acids/the_twenty.txt'
aacontent = open(aafile).read().splitlines()
aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)


def add_species():
    '''
    add species to the dataframe.
    :return:
    '''
    infile = outpath + '/' 'respairs_segment_notationx_len_merged.csv'
    pdbdir = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin/'
    df = pd.read_csv(infile)
    pdbids = df.pdbid.unique()
    lspecs = []
    hspecs = []
    aspecs = []
    amols = []
    notcons = []
    abchains = ['L', 'H']
    for i, row in df.iterrows():
        pdbid = row.pdbid
        pdbfile = pdbdir + pdbid + '.pdb'
        contents = open(pdbfile).read().splitlines()
        specieses = [item for item in contents if 'SPECIES' in item]
        if len(specieses)==3:
            lspec = specieses[0][26:].strip()
            hspec = specieses[1][26:].strip()
            aspec = specieses[2][26:].strip()
        else:
            chain_spec = dict((item[22].strip(), item[26:].strip()) for item in specieses)
            chain_labs = [item for item in contents if 'CHAIN' in item]
            chain_labs=  [item.split()[-2] for item in chain_labs if 'TYPE' not in item][:3]
            spe = [item.split(':')[1].strip() for item in specieses]
            for chain in chain_labs:
                if chain in chain_spec:
                    if chain == 'L':
                        lspec = chain_spec[chain]
                    elif chain == 'H':
                        hspec = chain_spec[chain]
                    else:
                        aspec = chain_spec[chain]
                else:
                    if chain == 'L':
                        lspec  = 'Not available'
                    elif chain == 'H':
                        hspec = 'Not available'
                    else:
                        aspec = 'Not available'
        lspecs.append(lspec.replace(',', ''))
        hspecs.append(hspec.replace(',', ''))
        aspecs.append(aspec.replace(',',''))
        if aspec == 'PHAGE #D': # sanity check
            print(pdbid)
    df['hspecies'] = hspecs
    df['lspecies'] = lspecs
    df['aspecies'] = aspecs
    # df['amolecule'] = amols
    outname = infile.split('.')[0] + '_species.csv'
    df.to_csv(outname, index=False)


def add_species_collapse():
    '''
    add species. merge species name
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species.csv'
    df = pd.read_csv(infile)
    aspecs = []
    for i, row in df.iterrows():
        aspec = str(row.aspecies)
        print(aspec)
        if 'CENTRUROIDES NOXIUS' in aspec:
            aspec = 'CENTRUROIDES NOXIUS'
        if 'DENGUE VIRUS' in aspec:
            aspec = 'DENGUE VIRUS'
        if 'ESCHERICHIA COLI' in aspec:
            aspec = 'ESCHERICHIA COLI'
        if 'INFLUENZA' in aspec:
            aspec = 'INFLUENZA'
        if 'HEPATITIS C VIRUS' in aspec:
            aspec = 'HEPATITIS C VIRUS'
        if 'HIV' in aspec or 'HUMAN IMMUNODEFICIENCY VIRUS' in aspec:
            aspec = 'HIV'
        if 'HERPESVIRUS' in aspec:
            aspec = 'HERPESVIRUS'
        if 'HUMAN RESPIRATORY SYNCYTIAL VIRUS' in aspec:
            aspec = 'HUMAN RESPIRATORY SYNCYTIAL VIRUS'
        if 'METHANOCALDOCOCCUS JANNASCHII' in aspec:
            aspec = 'METHANOCALDOCOCCUS JANNASCHII'
        if 'MUS MUSCULUS' in aspec:
            aspec = 'MUS MUSCULUS'
        if 'NEISSERIA MENINGITIDIS' in aspec:
            aspec = 'NEISSERIA MENINGITIDIS'
        if 'PLASMODIUM FALCIPARUM' in aspec:
            aspec = 'PLASMODIUM FALCIPARUM'
        if 'PLASMODIUM VIVAX' in aspec:
            aspec = 'PLASMODIUM VIVAX'
        if 'REOVIRUS' in aspec:
            aspec = 'REOVIRUS'
        if 'SACCHAROMYCES CEREVISIAE' in aspec:
            aspec = 'SACCHAROMYCES CEREVISIAE'
        if 'STAPHYLOCOCCUS AUREUS' in aspec:
            aspec = 'STAPHYLOCOCCUS AUREUS'
        if 'EBOLA' in aspec:
            aspec = 'EBOLA VIRUS'
        if 'CORONA' in aspec:
            aspec = 'CORONA VIRUS'
        aspecs.append(aspec)
    df['ag_species2'] = aspecs
    df = df.rename(columns = {'aspecies': 'ag_species'})
    outname = infile.split('.')[0] + '2.csv'
    df.to_csv(outname, index=False)

def add_full_sequence():
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
        # print(chain_list)
        # print(chain_dict[chain_list[0]])
        lseqs.append(chain_dict['L'])
        hseqs.append(chain_dict['H'])

    outdf = pd.DataFrame()
    outdf['pdbid'] = df.pdbid
    outdf['region'] = df.segment
    outdf['l_sequence'] = lseqs
    outdf['h_sequence'] = hseqs
    outdf['paratope'] = df.paratope
    print(outdf.head())
    outname = 'abdb_outfiles_2019/' + 'heavy_light_aaseq.csv'
    outdf.to_csv(outname , index=False)

def cdr_fr_motif_overlap():
    '''
    get motif overlap among segment
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
    df = pd.read_csv(infile)
    abchains = df.abchain.unique()
    segments = sorted(df.segment.unique().tolist())
    frs = sorted([item for item in segments if 'FR' in item])
    cdrs = sorted([item for item in segments if 'CDR' in item])
    motif_types = ['ab_motif', 'ag_motif']
    for motif_type in motif_types:
        data = []
        for i,segment in enumerate(segments):
            sdf = df[df.segment == segment]
            segment_abmotifs = getattr(sdf, motif_type).unique().tolist()
            coldata = [segment]
            for i2, segment2 in enumerate(segments):
                sdf2 = df[df.segment == segment2]
                segment_abmotifs2 = getattr(sdf2, motif_type).unique().tolist()
                overlap = set(segment_abmotifs) & set(segment_abmotifs2)
                lens = [len(segment_abmotifs), len(segment_abmotifs2)]
                max_len = max(lens)
                min_len = min(lens)
                # percent_overlap = round(len(overlap)/len(set(cdr_abmotifs +fr_abmotifs)),2)
                # percent_overlap = round(len(overlap)/len(fr_abmotifs),2)
                if i2 >= i:
                    percent_overlap = round(len(overlap)/max_len,2)
                else:
                    percent_overlap = round(len(overlap)/min_len, 2)
                percent_overlap = int(percent_overlap * 100)
                print(percent_overlap)
                coldata.append(percent_overlap)
            data.append(coldata)
        colnames = ['rownames'] + segments
        outdf  = pd.DataFrame(data = data, columns= colnames)
        outname = infile.split('.')[0] + '_%s_overlap_min_max.csv' % motif_type
        print(outname)
        outdf.to_csv(outname, index=False)


def get_abseq_full():
    '''
    Gets anibody sequence along with residue positions
    input:
    ------
    infiles = path to input files
    '''
    # filter for high res structure by matching with cutt of file
    cutoff_file = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments.csv'
    pdbdir = '../datasets/NR_LH_Protein_Martin/'
    pdbids = pd.read_csv(cutoff_file).pdbid.unique()
    infiles = [pdbdir+pdbid+'.pdb' for pdbid in pdbids]
    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    abchain = ['L', 'H']
    data = []
    lsegments = [('LFR1',(1,23)), ('CDR-L1', (24,34)), ('LFR2',(35,49)), ('CDR-L2',(50,56)),
                 ('LFR3',(57,88)),('CDR-L3',(89,97)), ('LFR4',(98,110))]
    lsegment_dict  = {}
    for lsegment in lsegments:
        newdict = dict([(i,lsegment[0]) for i in range(lsegment[1][0], lsegment[1][1]+1)])
        lsegment_dict.update(newdict)
    print(lsegment_dict)

    hsegments = [('HFR1',(1,30)), ('CDR-H1',(31,35)), ('HFR2',(36,49)), ('CDR-H2',(50,65)), ('HFR3',(66,94)),
                 ('CDR-H3',(95,102)), ('HFR4',(103,113))]
    hsegment_dict  = {}
    for hsegment in hsegments:
        newdict = dict([(i,hsegment[0]) for i in range(hsegment[1][0], hsegment[1][1]+1)])
        hsegment_dict.update(newdict)
    abchains = ['L', 'H']
    segment_dicts = dict([('L',lsegment_dict), ('H',hsegment_dict)])
    for fp in infiles:
        parser = PDB.PDBParser(PERMISSIVE=1)
        structure_name = fp.split('/')[-1].split('.')[0]
        structure = parser.get_structure(structure_name, fp)
        chains = structure.get_chains()
        for chain in chains:
            if chain.id in abchains:
                residues_list = chain.get_list()
                #sequence = ''.join([aadict[res.resname] for res in residues_list])
                #residue_number = '-'.join([''.join(str(item) for item in res.get_id()[1:]).strip() for res in residues_list])
                for res in residues_list:
                    het = res.get_id()[0] # residue type
                    if het == ' ': # take only standard residues
                        aa = res.resname
                        aa_single = aadict[aa]
                        resnum = res.get_id()[1]
                        resnumi = ''.join([str(item) for item in res.get_id()[1:]]).strip()
                        segment = segment_dicts[chain.id][resnum]
                        datum = (structure_name, aa, aa_single, resnum, resnumi, segment)
                        data.append(datum)
                else:
                    pass
    columns = ['pdbid', 'aa', 'aa_single', 'resnum', 'resnumi', 'segment']
    df = pd.DataFrame(data, columns=columns)
    df.to_csv('abdb_outfiles_2019/abdb_segment_absequence_full.csv',index=False)


def get_full_segment_vgene():
    '''
    gets full sequence instead of just intearcting residues plus the vgene section.
    :param infile:
    :return:
    '''
    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    infile = 'abdb_outfiles_2019/abdb_segment_absequence_full.csv'
    df = pd.read_csv(infile).iloc[:]
    data = []
    counter = 0
    print(df.head())
    for pdbid in df.pdbid.unique():
        pdbdf = df[df.pdbid==pdbid]
        segments = pdbdf.segment
        lvseq = '' # v gene is assumed FR1-FR3
        hvseq = '' # v gene is assumed FR1-FR3
        nonvs  = ['CDR-H3', 'CDR-L3', 'HFR4', 'LFR4']
        seqs = []
        pdbfile = '../datasets/NR_LH_Protein_Martin/%s.pdb' % pdbid
        pdbcontent = [item for item in open(pdbfile).readlines() if 'MOLECULE' in item or 'SPECIES' in item]
        pdbcontent = [item.split(':')[-2:] for item in pdbcontent]
        contents = open(pdbfile).read().splitlines()
        specieses = [item for item in contents if 'SPECIES' in item]
        molecules = [item for item in contents if 'MOLECULE' in item]
        if len(specieses)==3:
            lspec = specieses[0][26:].strip()
            hspec = specieses[1][26:].strip()
            aspec = specieses[2][26:].strip()
            lmol = molecules[0][26:].strip()
            hmol = molecules[1][26:].strip()
            amol = molecules[2][26:].strip()
        else:
            chain_spec = dict((item[22].strip(), item[26:].strip()) for item in specieses)
            chain_mol = dict((item[22].strip(), item[26:].strip()) for item in molecules)
            chain_labs = [item for item in contents if 'CHAIN' in item]
            chain_labs=  [item.split()[-2] for item in chain_labs if 'TYPE' not in item][:3]
            for chain in chain_labs:
                if chain in chain_spec:
                    if chain == 'L':
                        lspec = chain_spec[chain]
                    elif chain == 'H':
                        hspec = chain_spec[chain]
                    else:
                        aspec = chain_spec[chain]
                else:
                    if chain == 'L':
                        lspec  = 'Not available'
                    elif chain == 'H':
                        hspec = 'Not available'
                    else:
                        aspec = 'Not available'
            for chain in chain_mol:
                if chain in chain_mol:
                    if chain == 'L':
                        lmol = chain_mol[chain]
                    elif chain == 'H':
                        hmol = chain_mol[chain]
                    else:
                        amol = chain_mol[chain]
                else:
                    if chain == 'L':
                        lmol  = 'Not available'
                    elif chain == 'H':
                        hmol = 'Not available'
                    else:
                        amol = 'Not available'
        lspec = lspec.replace(',', '')
        hspec = hspec.replace(',', '')
        aspec = aspec.replace(',','')
        lmol = lmol.replace(',', '')
        hmol = hmol.replace(',', '')
        amol = amol.replace(',','')
        for segment in segments.unique():
            segdf = pdbdf[pdbdf.segment==segment]
            if segment not in nonvs and 'L' in segment:
                seq = ''.join(segdf.aa_single)
                lvseq += seq
            elif segment not in nonvs and 'H' in segment:
                seq = ''.join(segdf.aa_single)
                hvseq += seq
        for segment2 in segments.unique():
            segdf2 = pdbdf[pdbdf.segment==segment2]
            seq2 = ''.join(segdf2.aa_single)
            if 'L' in segment2:
                datum = [pdbid, segment2, seq2, lvseq, lmol, lspec,hmol,hspec,amol,aspec]
            elif 'H' in segment2:
                datum = [pdbid, segment2, seq2, hvseq, lmol, lspec,hmol,hspec,amol, aspec]
            data.append(datum)
    columns = ['pdbid', 'segment', 'segment_seq', 'vgene', 'lmolecule', 'lspecies', 'hmolecule', 'hspecies',
               'amolecule', 'aspecies']
    df = pd.DataFrame(data, columns=columns)
    outfile = infile[:-4] + '_vgene.csv'
    df.to_csv(outfile, index=False)
    print(outfile)
    print(df[df.segment=='CDR-H3'].segment_seq.value_counts())
    print(df[df.segment=='CDR-H3'].vgene.value_counts())



def annotate_vgenes():
    '''
    annoatate the pdb structures with the nearest vgenes (human and mouse only)
    :return:
    '''
    infile = 'abdb_outfiles_2019/abdb_segment_absequence_full_vgene.csv'
    vgene_file = 'abdb_outfiles_2019/imgt_vgenes.csv'
    df = pd.read_csv(infile)
    vgenedf = pd.read_csv(vgene_file)
    available_species = ['HOMO SAPIENS', 'MUS MUSCULUS']
    print(df.shape)
    # filter non human and non mouse
    df = df[(df.lspecies == 'HOMO SAPIENS') | (df.lspecies == 'MUS MUSCULUS')]
    df = df[(df.hspecies == 'HOMO SAPIENS') | (df.hspecies == 'MUS MUSCULUS')]
    pdbchains = []
    chains = []
    imgt_genenames = []
    imgt_species = []
    for i, row in df.iterrows():
        segment = row.segment
        if '-' in segment:
            pdbchain = segment[-2]
        else:
            pdbchain = segment[0]
        chain = pdbchain
        species_type = '%sspecies' % chain.lower()
        species = getattr(row, species_type)
        if species == 'HOMO SAPIENS':
            species = 'Homo sapiens' #IMGT formating
        elif species == 'MUS MUSCULUS':
            species = 'Mus musculus_C57BL/6' #IMGT formating
        if species == 'Mus musculus_C57BL/6' and chain  == 'L':
            chain = 'K'
        vdf = vgenedf[(vgenedf.chain == chain) & (vgenedf.species == species) & (vgenedf.functionality == 'F')]
        ref_vseqs = vdf.aa_seq
        vseq = row.vgene
        distances = []
        for i,ref_vseq in enumerate(ref_vseqs):
            ld = jellyfish.levenshtein_distance(vseq, ref_vseq)
            distances.append((ld,i))
        if pdbchain == 'CDR-H3' and chain =='K':
            print(row)
            sys.exit()
        print(pdbchain, chain, species)
        nearest_ref_index = sorted(distances)[0][1]
        nearest_ref_seq = ref_vseqs.iloc[nearest_ref_index]
        nearest_vgene = vdf.iloc[nearest_ref_index,]
        print(nearest_vgene)
        print(nearest_vgene.aa_seq)
        print(vseq)
        nearest_imgt_genename = nearest_vgene.imgt_genename
        print(nearest_imgt_genename)
        pdbchains.append(pdbchain)
        chains.append(chain)
        imgt_genenames.append(nearest_imgt_genename)
        imgt_species.append(species)
    df['pdbchain'] = pdbchains
    df['chains'] =  chains
    df['imgt_vgenename'] = imgt_genenames
    df['imgt_species'] = imgt_species
    print(df.pdbid.unique().shape)
    outfile = infile[:-4] + '_imgt_vgene.csv'
    print(outfile)
    df.to_csv(outfile, index=False)

def merge_notationx_vgene():
    '''
    ad motif to vgene file
    :return:
    '''
    nfile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
    vfile = 'abdb_outfiles_2019/abdb_segment_absequence_full_vgene_imgt_vgene.csv'
    vdf = pd.read_csv(vfile)
    ndf = pd.read_csv(nfile)
    vgenes = []
    rows = []
    imgt_speciess = []
    available_pdbids = vdf.pdbid.unique()
    for i, row in ndf.iterrows():
        pdbid = row.pdbid
        if pdbid in available_pdbids:
            segment = row.segment
            vndf = vdf[(vdf.pdbid == pdbid) & (vdf.segment == segment)]
            vgene = vndf.iloc[0].imgt_vgenename
            vgenes.append(vgene)
            rows.append(row)
            imgt_species = vndf.iloc[0].imgt_species
            imgt_speciess.append(imgt_species)
    outdf = pd.DataFrame(rows)
    outdf['imgt_vgenename'] = vgenes
    outdf['imgt_species'] = imgt_speciess
    print(outdf.pdbid.unique().shape)
    outname = nfile.split('.')[0] + '_imgt_vgene.csv'
    outdf.to_csv(outname, index=False)


# run stuff
# add_species()
# add_species_collapse()
# add_full_sequence()
# cdr_fr_motif_overlap()
# get_abseq_full()
# get_full_segment_vgene()
# annotate_vgenes()
# merge_notationx_vgene()



