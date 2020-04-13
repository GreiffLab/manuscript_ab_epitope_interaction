# import stuff
import sys
import pandas as pd
import math

#set display
pd.set_option('display.max_column', None)

def get_residue_logodd():
    '''
    get residue pairs log odds. follows kunik & ofran (https://academic.oup.com/peds/article/26/10/599/1510690)
    annotate by regions and motifs...
    :return:
    '''
    infile =  'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments.csv'
    df = pd.read_csv(infile)
    print(df.head())
    segments = df.segment.unique()
    resnames = df.abres.unique()
    resname_pairs = [['%s_%s' % (res1, res2) for res1 in resnames] for res2 in resnames]
    resname_pairs = sum(resname_pairs,[])
    print(resnames)
    print(resname_pairs)
    print(len(resname_pairs))
    print(segments)
    data = []
    for segment in segments[:]:
        sdf =  df[df.segment == segment]
        print(sdf.shape, segment, df.shape)
        for res1 in resnames:
            for res2 in resnames:
                nres1 = sdf[sdf.abres==res1].shape[0]
                nres2 = sdf[sdf.agres==res2].shape[0]
                print(nres1, nres2)
                pdf = sdf[(sdf.abres==res1) & (sdf.agres==res2)]
                npair  = pdf.shape[0]
                print(pdf)
                nsegment = float(sdf.shape[0])
                ppair = npair/nsegment
                pres1 = nres1/nsegment
                pres2 = nres2/nsegment
                if pres1*pres2 > 0:
                    odd = ppair/(pres1*pres2)
                else:
                    strodd = 'NA'
                print(odd, 'odd')
                if odd > 0:
                    lodd = math.log2(odd)
                elif odd == 0 or strodd == 'NA':
                    lodd = 'NA'
                print(lodd,'lodd')
                datum = [res1, res2, segment, odd, lodd]
                data.append(datum)
    outname = infile.split('.')[0] + '_residue_contact_odds.csv'
    colnames  = ['abres', 'agres', 'region', 'odd', 'logodd']
    outdf = pd.DataFrame(data, columns=colnames)
    print(outdf.shape[0]/len(segments))
    print(outname)
    outdf.to_csv(outname, index=False)



def get_residue_logodd_ppi():
    '''
    get residue pairs log odds. follows kunik & ofran (https://academic.oup.com/peds/article/26/10/599/1510690)
    annotate by regions and motifs...
    :return:
    '''
    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    print(aadict)
    aadictr = dict([(value, key) for key, value in aadict.items()])
    infile =  'abdb_outfiles_2019/threedid_no_iglike.csv'
    df = pd.read_csv(infile)
    infile2 = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300.csv'
    df2 = pd.read_csv(infile2)
    dff = df[df.pdbid.isin(df2.pdbid)]
    print(dff.head())
    intertypes = dff.intertype.unique()
    print(intertypes)
    resnames = df.res1.unique()
    print(resnames)
    # resnames = resnames[:2]
    data = []
    for intertype in intertypes:
        sdf = dff[dff.intertype == intertype]
        print(sdf.shape)
        nsegment = float(sdf.shape[0])
        for res1 in resnames:
            for res2 in resnames:
                nres1 = sdf[sdf.res1==res1].shape[0]
                nres2 = sdf[sdf.res2==res2].shape[0]
                print(nres1, nres2)
                pdf = sdf[(sdf.res1==res1) & (sdf.res2==res2)]
                npair  = pdf.shape[0]
                print(pdf.shape, '# of pair')
                ppair = npair/nsegment
                pres1 = nres1/nsegment
                pres2 = nres2/nsegment
                if pres1*pres2 > 0:
                    odd = ppair/(pres1*pres2)
                else:
                    strodd = 'NA'
                print(odd, 'odd')
                if odd > 0:
                    lodd = math.log2(odd)
                elif odd == 0 or strodd == 'NA':
                    lodd = 'NA'
                print(lodd,'lodd')
                res1three = aadictr[res1]
                res2three = aadictr[res2]
                datum = [res1three, res2three, intertype, odd, lodd]
                data.append(datum)
    outname = infile2.split('.')[0] + '_residue_contact_odds.csv'
    colnames  = ['abres', 'agres', 'region', 'odd', 'logodd']
    outdf = pd.DataFrame(data, columns=colnames)
    print(outdf.shape[0]/len(intertypes))
    print(outname)
    outdf.to_csv(outname, index=False)

## run stuff
# get_residue_logodd()
get_residue_logodd_ppi()