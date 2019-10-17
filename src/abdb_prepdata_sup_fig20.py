# import stuff
import sys
import pandas as pd
import jellyfish
import os


pd.set_option('display.max_column', None)

def gap_in_seq():
    '''
    insert gaps in within paratope/epitope sequences
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber.csv'
    df = pd.read_csv(infile)
    df = df.dropna(subset=['paratope', 'epitope'])
    print(df.head())
    abgapmotifs = []
    abgapmotifs2 = []
    abgapmotifs3 = []
    aggapmotifs = []
    aggapmotifs2 = []
    aggapmotifs3 = []
    for i, row in df.iterrows():
        abgaps = row.ab_motif.split('X')
        print(abgaps, row.ab_motif)
        aggaps = row.ag_motif.split('X')
        print(row.paratope)
        pchars = list(row.paratope)
        echars = list(row.epitope)
        abgapmotif  = ''
        abgapmotif2 = ''
        abgapmotif3 = ''
        for res,gap in zip(pchars, abgaps[1:]):
            abgapmotif += res + gap
            if gap != '':
                abgapmotif2 += res + '-'*int(gap)
                abgapmotif3 += res + '-'
            else:
                abgapmotif2 += res + gap
                abgapmotif3 += res + gap
        aggapmotif  = ''
        aggapmotif2 = ''
        aggapmotif3 = ''
        for res,gap in zip(echars, aggaps[1:]):
            aggapmotif += res + gap
            if gap != '':
                aggapmotif2 += res + '-'*int(gap)
                aggapmotif3 += res + '-'
            else:
                aggapmotif2 += res + gap
                aggapmotif3 += res + gap
        print(aggapmotif, aggapmotif2,row.ag_motif, row.epitope)
        aggapmotifs.append(aggapmotif)
        aggapmotifs2.append(aggapmotif2)
        aggapmotifs3.append(aggapmotif3)
        abgapmotifs.append(abgapmotif)
        abgapmotifs2.append(abgapmotif2)
        abgapmotifs3.append(abgapmotif3)
    df['aggapmotif'] = aggapmotifs
    df['aggapmotif2'] = aggapmotifs2
    df['aggapmotif3'] = aggapmotifs3
    df['abgapmotif'] = abgapmotifs
    df['abgapmotif2'] = abgapmotifs2
    df['abgapmotif3'] = abgapmotifs3
    outname = infile.split('.')[0] + '_phil.csv'
    df.to_csv(outname, index=False)

def prep_data_branchld(param, param2):
    '''
    preps data for branch ld
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil_pc.csv'
    df = pd.read_csv(infile)
    df = df[(df.plen > 1) & (df.epitope_len >1)]
    # df = df[(df.epitope_len > 1)]
    branches = getattr(df, param).unique()
    print(len(branches))
    sizes = []
    avelds = []
    seens = []
    for i, branch in enumerate(branches):
        rdf = df[getattr(df,param) == branch]
        epitopes = getattr(rdf, param2)
        # if len(epitopes) == 1:
        for epi in epitopes:
            lds = []
            for epi2 in epitopes:
                # print(epi, epi2, 'hey')
                ld = jellyfish.levenshtein_distance(epi, epi2)/max(len(epi), len(epi2))
                if ld != 0:
                    ld = 1
                lds.append(ld)
            # print(lds)
            aveld = sum(lds)/len(lds)
            avelds.append(aveld)
            # print(aveld)
            # print(df.shape)
    avelds = sum(avelds)/len(avelds)
    print(avelds)

def physchem_encoding():
    '''
    encode paratope motif with physchem props
    :return:
    grouping in accordance to the R package Peptides and Rice,  Peter,  Ian Longden,  and Alan Bleasby.
    "EMBOSS: the European molecular biology opensoftware suite." Trends in genetics 16.6 (2000): 276-277.
    Polar (P), Aromatic (R), Non-polar (N), Charged (C)
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil.csv'
    df = pd.read_csv(infile)
    print(df.head())
    aapcdict = {
        'D':'C',
        'E':'C',
        'H':'C',
        'R':'C',
        'K':'C',
        'S':'P',
        'T':'P',
        'N':'P',
        'Q':'P',
        'G':'N',
        'A':'N',
        'V':'N',
        'L':'N',
        'I':'N',
        'P':'N',
        'M':'N',
        'C':'N',
        'F':'R',
        'W':'R',
        'Y':'R',
        '-':'-'
    }
    parpcmotifs = []
    epipcmotifs = []
    for i, row in df.iterrows():
        parmotif = row.abgapmotif3
        epimotif = row.aggapmotif3
        print(parmotif)
        print(epimotif)
        parpcmotif = ''.join([aapcdict[res] for res in parmotif])
        epipcmotif = ''.join([aapcdict[res] for res in epimotif])
        print(parpcmotif, epipcmotif)
        parpcmotifs.append(parpcmotif)
        epipcmotifs.append(epipcmotif)
    df['parpcmotif'] = parpcmotifs
    df['epipcmotif'] = epipcmotifs
    print(df)
    outfile = infile.split('.')[0] + '_pc.csv'
    print(outfile)
    df.to_csv(outfile, index=False)

def resgapmotif_dataset(datasetdir, para, epi):
    '''
    prepare dl dataset for residue-gap encodeing (aggregate encoding)
    drop single chars
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil_pc.csv'
    df =  pd.read_csv(infile)
    print(df.shape)
    # df = df[(df.plen > 1) & (df.epitope_len >1)]
    print(df.shape)
    # print(df.head())
    # sys.exit()
    datasetdir = '/Users/rahmadakbar/greifflab/aims/aimugen/dl/dataset_%s' % datasetdir
    os.system('mkdir %s' % datasetdir)
    outname = '%s/paraepi.tsv' % datasetdir
    outname2 = '%s/epipara.tsv' % datasetdir
    paras = getattr(df, para)
    epis = getattr(df, epi)
    outcontent = '\n'.join(['\t'.join(item) for item in zip(paras, epis)])
    outcontent2 = '\n'.join(['\t'.join(item) for item in zip(epis, paras)])
    outfile  = open(outname, 'w')
    outfile.write(outcontent)
    outfile2 = open(outname2, 'w')
    outfile2.write(outcontent2)


#run stuff
# gap_in_seq()
# prep_data_branchld('abgapmotif', 'aggapmotif')
# prep_data_branchld('abgapmotif3', 'aggapmotif3')
# prep_data_branchld('ab_motif', 'ag_motif')
# prep_data_branchld('paratope', 'epitope')
# physchem_encoding()
# prep_data_branchld('parpcmotif', 'epipcmotif')
resgapmotif_dataset('ressingle', 'abgapmotif3', 'aggapmotif3')
# resgapmotif_dataset('pc', 'parpcmotif', 'epipcmotif')
