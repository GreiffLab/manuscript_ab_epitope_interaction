# import stuff
import sys
import pandas as pd
import jellyfish
import os
from find_files import find_files as fifi

pd.set_option('display.max_column', None)


def sl_dl_summary():
    '''
    merged summary stats from dl and sl models
    :return:
    '''
    infiles = fifi('../sl/results', 'summary')
    infiles2 = fifi('../sl/results', 'summary')
    print(len(infiles))
    # exclude branch files
    infiles = [item for item in infiles if 'ppi/' in item]
    infiles = infiles  + [item for item in infiles2 if 'res/' in item]
    print(len(infiles))
    print(infiles)
    # sys.exit()
    exacts = []
    norms = []
    for infile in infiles:
        df = pd.read_csv(infile)
        # print(df)
        print(infile)
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
        elif '_res' in infile:
            data_tag = 'agg'
            if 'epitope' in infile:
                use_case_tag = 'agg_epipara'
            else:
                use_case_tag = 'agg_paraepi'
        else:
            data_tag = 'seq'
            if 'epitope' in infile:
                use_case_tag = 'seq_epipara'
            else:
                use_case_tag = 'seq_paraepi'
        if 'ppi' in infile:
            source_tag = 'ppi'
            exp_tag = 'exp'
            if 'randomized' in infile:
                exp_tag = 'control'
        elif '_res' in infile:
            source_tag = 'abdb'
            exp_tag = 'exp'
            if 'randomized' in infile:
                exp_tag = 'control'
        else:
            source_tag = 'abdb'
            exp_tag = 'exp'
            if 'randomized' in infile:
                exp_tag = 'control'
        if 'exact' in infile:
            base_types = ['marginal_proba', 'cond_proba', 'cond_proba_with_prior']
            for base_type in base_types:
                exp_tag_bt = exp_tag + '_' + base_type
                print(exp_tag_bt)
                ld = df[df.approach == base_type]
                print(ld)
                error = ld.error_mean.iloc[0]
                error_sd = ld.error_standard_deviation.iloc[0]
                print(error, error_sd,'hey')
                exact = [use_case_tag, data_tag, exp_tag_bt,error, error_sd, 10, source_tag]
                exacts.append(exact)
                print(exact)
        elif 'LD' in infile:
            base_types = ['marginal_proba', 'cond_proba', 'cond_proba_with_prior']
            for base_type in base_types:
                exp_tag_bt = exp_tag + '_' + base_type
                print(exp_tag_bt)
                ld = df[df.approach == base_type]
                error = ld.error_mean.iloc[0]
                error_sd = ld.error_standard_deviation.iloc[0]
                print(error, error_sd,'hey')
                norm = [use_case_tag, data_tag, exp_tag_bt,error, error_sd, 10, source_tag]
                norms.append(norm)
                print(norm)
    print(len(norms), len(exacts))
    colnames = ['use_case','data_tag','exp_tag','repldnormmea','repldnormse',
                'ldnormreps','repldexactmea','repldexactse','ldexactreps', 'source']
    outdata = [item1 + item2[-3:] for item1, item2 in zip(norms, exacts)]
    outdf = pd.DataFrame(outdata, columns=colnames)
    print(outdf.shape)
    # sys.exit()
    evalsumfile = 'abdb_outfiles_2019/eval_summary.csv'
    dfeval = pd.read_csv(evalsumfile)
    print(dfeval)
    # sys.exit()
    aggevalsumfile = 'abdb_outfiles_2019/agg_eval_summary.csv'
    aggdfeval = pd.read_csv(aggevalsumfile)
    print(aggdfeval)
    # sys.exit()
    mergeddf = pd.concat([aggdfeval,dfeval, outdf])
    print(mergeddf.head())
    print(mergeddf.exp_tag)
    # sys.exit()
    outname = 'abdb_outfiles_2019/sl_dl_evalsummary_ppi.csv'
    cats = []
    # exp_tag2s = []
    # exp_tag_dict = {'control': 'control',
    #                 'exp': 'Imm deep NMT',
    #                 'exp_base_majority': 'Imm shallow majority',
    #                 'exp_base_mapping': 'Imm shallow mapping',
    #                 'exp_base_proba': 'Imm shallow proba',
    #                 'exp_base_ppi_majority': 'Non-imm shallow majority',
    #                 'exp_base_ppi_mapping': 'Non-imm shallow mapping',
    #                 'exp_base_ppi_proba': 'Non-imm shallow proba'}
    for i,row in mergeddf.iterrows():
        cat = row.exp_tag.split('_')[0]
        cats.append(cat)
    mergeddf['category'] = cats
    print(mergeddf)
    # sys.exit()
    mergeddf.to_csv(outname, index=False)


#run stuff

sl_dl_summary()
