# import stuff
# for foldx, do run_foldx.py

import pandas as pd
import os
import sys


def renumber_foldx_outresidues():
    '''
    renumber foldx residues in accordance with abdb-martin
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil_pc.csv'
    df = pd.read_csv(infile)
    pdbids = df.pdbid.unique()
    n = len(pdbids)
    print(n)
    nrpath = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin'
    foldxpath = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin_foldx'
    lsegments = [('LFR1', (1, 23)), ('CDR-L1', (24, 34)), ('LFR2', (35, 49)), ('CDR-L2', (50, 56)),
                 ('LFR3', (57, 88)), ('CDR-L3', (89, 97)), ('LFR4', (98, 110))]
    lsegment_dict = {}
    for lsegment in lsegments:
        newdict = dict([(i, lsegment[0]) for i in range(lsegment[1][0], lsegment[1][1] + 1)])
        lsegment_dict.update(newdict)
    print(lsegment_dict)

    hsegments = [('HFR1', (1, 30)), ('CDR-H1', (31, 35)), ('HFR2', (36, 49)), ('CDR-H2', (50, 65)), ('HFR3', (66, 94)),
                 ('CDR-H3', (95, 102)), ('HFR4', (103, 113))]
    hsegment_dict = {}
    for hsegment in hsegments:
        newdict = dict([(i, hsegment[0]) for i in range(hsegment[1][0], hsegment[1][1] + 1)])
        hsegment_dict.update(newdict)
    print(lsegment_dict)
    print(hsegment_dict)
    data = []
    for pdbid in pdbids:
        nfile = nrpath + '/' + pdbid + '.pdb'
        ffile  = foldxpath + '/' + pdbid + '_AS.fxout'
        fcontents = open(ffile).read().splitlines()
        tracker = 0
        chain_tracker = -1
        chains = ['L', 'H', 'A']
        for fcontent in fcontents:
            parts = fcontent.split()
            resnum = int(parts[1])
            print(parts)
            dg = abs(float(parts[-1]))
            if resnum <= 1:
                if chain_tracker == -1:
                    chain_tracker += 1
                    chain = chains[chain_tracker]
                elif chain_tracker==0:
                    chain_tracker = 1
                    chain = chains[chain_tracker]
                    print(chain)
                    # sys.exit()
            if chain == 'L' and 1<=resnum <=110:
                res  = parts[0]
                segment = lsegment_dict[resnum]
                datum = [pdbid, res, resnum, chain, segment, dg]
                print(datum)
            elif chain == 'H' and 1<=resnum <=113:
                res  = parts[0]
                segment = hsegment_dict[resnum]
                datum = [pdbid, res, resnum, chain, segment, dg]
                print(datum)
            data.append(datum)
    colnames = ['pdbid', 'abres', 'abresnum', 'abchain', 'segment', 'deltag']
    df = pd.DataFrame(data, columns=colnames)
    df = df[df.deltag > 10]
    print(df.head())
    print(df.shape)
    outfile = 'abdb_outfiles_2019/foldx_deltag.csv'
    df.to_csv(outfile, index=False)


def renumber_foldx_splitoutresidues():
    '''
    renumber foldx residues in accordance with abdb-martin
    :return:
    '''
    foldxpath = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin_foldx/chain_split'
    outfiles = os.listdir(foldxpath)
    print(outfiles)
    outfiles = [os.path.join(foldxpath,item) for item in outfiles if '.fxout' in item]
    lsegments = [('LFR1', (1, 23)), ('CDR-L1', (24, 34)), ('LFR2', (35, 49)), ('CDR-L2', (50, 56)),
                 ('LFR3', (57, 88)), ('CDR-L3', (89, 97)), ('LFR4', (98, 110))]
    lsegment_dict = {}
    for lsegment in lsegments:
        newdict = dict([(i, lsegment[0]) for i in range(lsegment[1][0], lsegment[1][1] + 1)])
        lsegment_dict.update(newdict)
    print(lsegment_dict)

    hsegments = [('HFR1', (1, 30)), ('CDR-H1', (31, 35)), ('HFR2', (36, 49)), ('CDR-H2', (50, 65)), ('HFR3', (66, 94)),
                 ('CDR-H3', (95, 102)), ('HFR4', (103, 113))]
    hsegment_dict = {}
    for hsegment in hsegments:
        newdict = dict([(i, hsegment[0]) for i in range(hsegment[1][0], hsegment[1][1] + 1)])
        hsegment_dict.update(newdict)
    print(lsegment_dict)
    print(hsegment_dict)
    data = []
    excepted  = 0
    for outfile in outfiles:
        if '_H_' in outfile:
            chain = 'H'
            pdbid = outfile.split('/')[-1].split('_')[0]
            fcontents = open(outfile).read().splitlines()
            maxn = 1
            seens = range(1, maxn)
            for i,content in enumerate(fcontents):
                parts = content.split()
                res = parts[0]
                resnum = int(parts[1])
                # deltag = abs(float(parts[-1]))
                deltag = float(parts[-1])
                if resnum not in seens and 0<resnum <= 113:
                    segment = hsegment_dict[resnum]
                    datum = [pdbid, res, resnum, chain, deltag, segment]
                    data.append(datum)
                    print(datum)
                    # seens.append(resnum)
                    maxn=resnum
        if '_L_' in outfile:
            chain = 'L'
            pdbid = outfile.split('/')[-1].split('_')[0]
            fcontents = open(outfile).read().splitlines()
            maxn = 1
            seens = range(1,maxn)
            print(outfile)
            for i,content in enumerate(fcontents):
                parts = content.split()
                res = parts[0]
                resnum = int(parts[1])
                # deltag = abs(float(parts[-1]))
                deltag = float(parts[-1])
                if resnum not in seens and 0 < resnum <=110:
                    segment = lsegment_dict[resnum]
                    datum = [pdbid, res, resnum, chain, deltag, segment]
                    data.append(datum)
                    print(datum)
                    print(resnum)
                    # seens.append(resnum)
                    maxn=resnum
    colnames=  ['pdbid', 'abres', 'abresnum', 'abchain', 'deltag', 'segment']
    outdf = pd.DataFrame(data, columns=colnames)
    print(outdf.head())
    outfile = 'abdb_outfiles_2019/foldx_deltag_split_chain.csv'
    outdf.to_csv(outfile, index=False)

def parse_beatoutfile():
    '''
    parses beatmusic outfiles for comparison to foldx
    :return:
    '''
    infile = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin_foldx/chain_split/beatmusic/1a14_1_h.pdb.beat.txt'
    contents = open(infile).read().splitlines()
    lsegments = [('LFR1', (1, 23)), ('CDR-L1', (24, 34)), ('LFR2', (35, 49)), ('CDR-L2', (50, 56)),
                 ('LFR3', (57, 88)), ('CDR-L3', (89, 97)), ('LFR4', (98, 110))]
    lsegment_dict = {}
    for lsegment in lsegments:
        newdict = dict([(i, lsegment[0]) for i in range(lsegment[1][0], lsegment[1][1] + 1)])
        lsegment_dict.update(newdict)
    print(lsegment_dict)

    hsegments = [('HFR1', (1, 30)), ('CDR-H1', (31, 35)), ('HFR2', (36, 49)), ('CDR-H2', (50, 65)), ('HFR3', (66, 94)),
                 ('CDR-H3', (95, 102)), ('HFR4', (103, 113))]
    hsegment_dict = {}
    for hsegment in hsegments:
        newdict = dict([(i, hsegment[0]) for i in range(hsegment[1][0], hsegment[1][1] + 1)])
        hsegment_dict.update(newdict)
    data = []
    for content in contents:
        parts = content.split()
        print(parts)
        pdbid = parts[0].split('_')[0]
        abchain = parts[2]
        try:
            abresnum = int(parts[3])
        except:
            abres = int(parts[3][:-2])
        abres = parts[4]
        mutres = parts[5]
        deltag = float(parts[7])
        if abchain == 'H':
            segment = hsegment_dict[abresnum]
        if mutres == 'A':
            datum = [pdbid, abres, abresnum, abchain, deltag, segment]
            data.append(datum)
    print(data)
    colnames=  ['pdbid', 'abres', 'abresnum', 'abchain', 'deltag', 'segment']
    outdf = pd.DataFrame(data, columns=colnames)
    print(outdf.head())
    outfile = 'abdb_outfiles_2019/beatmusic_deltag_split_chain.csv'
    outdf.to_csv(outfile, index=False)


# runs stuff
# renumber_foldx_outresidues()
# renumber_foldx_splitoutresidues()
parse_beatoutfile()
