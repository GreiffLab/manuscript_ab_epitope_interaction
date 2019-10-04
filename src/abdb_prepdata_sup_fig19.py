# import stuf

import sys
import pandas as pd
import os
from find_files import find_files as fifi
# from abdb_prepdata_main_fig6 import make_tabsep_ppi
pd.set_option('display.max_column', None)

def branch_dataset():
    '''
    output sequence dataset for largest paratope branches (XXX)
    :return:
    '''
    nodesfile = 'abdb_outfiles_2019/paratope_epitope_internet_nodes.csv'
    df = pd.read_csv(nodesfile)
    n = 3
    topmotifs = df.id.iloc[:n]
    abagfile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber.csv'
    abagdf = pd.read_csv(abagfile)
    print(abagdf.info())
    for topmotif in topmotifs:
        tdf = abagdf[abagdf.ab_motif == topmotif]
        print(tdf.info())
        datasetdir = '/Users/rahmadakbar/greifflab/aims/aimugen/dl/dataset_%s' % topmotif
        os.system('mkdir %s' % datasetdir)
        outname = '%s/paraepi.tsv' % datasetdir
        outname2 = '%s/epipara.tsv' % datasetdir
        paras = tdf.paratope
        epis = tdf.epitope
        outcontent = '\n'.join(['\t'.join(item) for item in zip(paras, epis)])
        outcontent2 = '\n'.join(['\t'.join(item) for item in zip(epis, paras)])
        outfile  = open(outname, 'w')
        outfile.write(outcontent)
        outfile2 = open(outname2, 'w')
        outfile2.write(outcontent2)



def uniquepair_ppi_dataset():
    '''
    get unique pairs from ppi dataset
    :return:
    '''
    infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv'
    df = pd.read_csv(infile)
    pairs = []
    for i,row in df.iterrows():
        seq1 = row.sequence1
        seq2 = row.sequence2
        pair = '-'.join([seq1,seq2])
        pairs.append(pair)
    print(pairs)
    df['seqpair'] = pairs
    print(df.shape)
    df = df.drop_duplicates(subset='seqpair')
    print(df.shape)
    outfile = infile.split('.')[0] + '_uniqueseqpair.csv'
    df.to_csv(outfile, index=False)
    outseq1 = '../dl/dataset_ppi_uniqueseqpair/paraepi.tsv'
    outseq2 = '../dl/dataset_ppi_uniqueseqpair/epipara.tsv'
    outcontent1 = '\n'.join(['\t'.join([seq1,seq2]) for seq1,seq2 in zip(df.sequence1, df.sequence2)])
    outcontent2 = '\n'.join(['\t'.join([seq2,seq1]) for seq1,seq2 in zip(df.sequence1, df.sequence2)])
    print(outcontent1)
    print(outcontent2)
    outfile1 = open(outseq1, 'w')
    outfile1.write(outcontent1)
    outfile2 = open(outseq2, 'w')
    outfile2.write(outcontent2)



def sl_dl_summary():
    '''
    merged summary stats from dl and sl models
    :return:
    '''
    infiles = fifi('../sl', 'summary')
    branch_infiles = [item for item in infiles if 'XXX' in item]
    exacts = []
    norms = []
    for infile in branch_infiles:
        df = pd.read_csv(infile)
        print(df)
        # sys.exit()
        print(infile)
        branch_tag = infile.split('/')[3].split('_')[1]
        print(branch_tag)
        base_types = ['marginal_proba', 'cond_proba', 'cond_proba_with_prior']
        if 'motif' in infile:
            data_tag = 'motif'
            if 'epitope' in infile and 'pos' in infile:
                use_case_tag = 'motif_epiparapos'
            elif 'epitope' in infile and 'pos' not in infile:
                use_case_tag = 'motif_epipara'
            elif 'paratope' in infile and 'pos' in infile:
                use_case_tag = 'motif_paraepipos'
            else:
                use_case_tag = 'motif_paraepi'
        else:
            data_tag = 'seq'
            if 'epitope' in infile:
                use_case_tag = 'seq_epipara'
            else:
                use_case_tag = 'seq_paraepi'
        if 'ppi' in infile:
            exp_tag = 'exp_base_ppi'
        else:
            exp_tag = 'exp_base'
        if 'randomized' in infile:
            control_tag = 'control'
            print(infile, control_tag)
        else:
            control_tag = 'exp'
        if 'exact' in infile:
            for base_type in base_types:
                exp_tag_bt = exp_tag + '_' + base_type + '_' + control_tag
                print(exp_tag_bt)
                ld = df[df.approach == base_type]
                print(ld)
                error = ld.error_mean.iloc[0]
                error_sd = ld.error_standard_deviation.iloc[0]
                print(error, error_sd,'hey')
                exact = [use_case_tag, data_tag, exp_tag_bt, branch_tag, error, error_sd, 10]
                exacts.append(exact)
                print(exact)
        elif 'LD' in infile:
            for base_type in base_types:
                exp_tag_bt = exp_tag + '_' + base_type
                exp_tag_bt = exp_tag + '_' + base_type + '_' + control_tag
                print(exp_tag_bt)
                ld = df[df.approach == base_type]
                error = ld.error_mean.iloc[0]
                error_sd = ld.error_standard_deviation.iloc[0]
                print(error, error_sd,'hey')
                norm = [use_case_tag, data_tag, exp_tag_bt,branch_tag,error, error_sd, 10]
                norms.append(norm)
                print(norm)
    print(len(norms), len(exacts))
    colnames = ['use_case','data_tag','exp_tag', 'branch', 'repldnormmea','repldnormse',
                'ldnormreps','repldexactmea','repldexactse','ldexactreps']
    outdata = [item1 + item2[-3:] for item1, item2 in zip(norms, exacts)]
    outdf = pd.DataFrame(outdata, columns=colnames)
    print(outdf)
    evalsumfile = 'abdb_outfiles_2019/branch_eval_summary.csv'
    dfeval = pd.read_csv(evalsumfile)
    print(dfeval.head())
    print(outdf.head())
    mergeddf = pd.concat([dfeval, outdf])
    print(mergeddf.tail())
    # sys.exit()
    outname = 'abdb_outfiles_2019/branch_sl_dl_evalsummary.csv'
    cats = []
    exp_tag2s = []
    exp_tag_dict = {'control': 'Imm deep NMT ctrl',
                    'exp': 'Imm deep NMT',
                    'exp_base_marginal_proba_exp': 'Imm shallow marg',
                    'exp_base_cond_proba_exp': 'Imm shallow cond',
                    'exp_base_cond_proba_with_prior_exp': 'Imm shallow condpr',
                    'exp_base_marginal_proba_control': 'Imm shallow marg ctrl',
                    'exp_base_cond_proba_control': 'Imm shallow cond ctrl',
                    'exp_base_cond_proba_with_prior_control': 'Imm shallow condpr ctrl',
                    }
    for i,row in mergeddf.iterrows():
        print(row)
        exp_tag2 = exp_tag_dict[row.exp_tag]
        print(exp_tag2)
        cat = exp_tag2.split()[0]
        cats.append(cat)
        exp_tag2s.append(exp_tag2)
    mergeddf['exp_tag2'] = exp_tag2s
    mergeddf['category'] = cats
    print(mergeddf)
    # sys.exit()
    mergeddf.to_csv(outname, index=False)



# run stuff
# branch_dataset()
# uniquepair_ppi_dataset()
sl_dl_summary()




