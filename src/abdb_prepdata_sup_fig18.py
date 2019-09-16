# import stuff
import pandas as pd
import sys
import random


def victor_cumulative_curve():
    '''
    prepdata for victor's cumulative curve
    :return:
    '''
    infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv'
    df = pd.read_csv(infile)
    print(df.info())
    data = []
    n = 1000
    pdbids = df.pdbid1.unique()
    print(len(pdbids))
    print(df.head())
    #skips some structure so we can see the plot
    skip = 30
    indices = range(0,len(pdbids), skip)
    abag = ['gap_pattern1', 'gap_pattern2']
    abag2 = ['ab', 'ag']
    for motiftype, source in zip(abag, abag2):
        # source = motiftype.split('_')[0]
        for run in range(n):
            random.shuffle(pdbids)
            current_motif = []
            for i,pdbid in enumerate(pdbids):
                if i in indices:
                    pdbdf = df[df.pdbid1 == pdbid]
                    strpdbidx = 's%s' % str(i+2)
                    pdbidx = i+2
                    pdbmotif = getattr(pdbdf,motiftype).tolist()
                    if i == 0:
                        pdbid2 = pdbids[1]
                        pdbdf2 = df[df.pdbid1 == pdbid2]
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
    infile_tag = infile.split('.')[0]
    outname = '%s_motif_coverage_skip%s.csv' % (infile_tag,skip)
    print(outname)
    covdf.to_csv(outname, index=False)



# run stuff
# victor_cumulative_curve()

