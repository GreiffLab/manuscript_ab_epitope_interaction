# import stuff
import sys
import pandas as pd
import jellyfish
import os
from mpi4py import MPI
import numpy as np

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


def get_ppi_resnum():
    '''
    get ppi resnum from the original thredid flat file
    :return:
    '''
    flatfile = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/3did/3did_flat.txt'
    flatcontents = open(flatfile).read().split('//')
    infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv'
    df = pd.read_csv(infile)
    print(df.head())
    resnum1s = []
    resnum2s = []
    inpdbids = df.pdbchainpair1.unique()
    counter = 0
    outdfs  = []
    for flatcontent in flatcontents[:]:
        contents = flatcontent.split('#=3D')[1:]
        for content in contents:
            parts = content.splitlines()[0].split()
            pdbid = parts[0].upper()
            chain1 = parts[1][0].upper()
            chain2 = parts[2][0].upper()
            # print(pdbid, chain1, chain2)
            inpdbid = pdbid + '_' + chain1+ chain2
            if inpdbid in inpdbids:
                rescontent = ['res1\tres2\tresnum1\tresnum2\tssormm'] + content.splitlines()[1:]
                icontent = '\n'.join(rescontent)
                iocontent = StringIO(icontent)
                dfr = pd.read_csv(iocontent, sep='\t', )
                dfr1 = dfr.drop_duplicates(subset='resnum1')
                # dfr2 = dfr.sort_values(by = 'resnum2').drop_duplicates(subset='resnum2')
                dfr2 = dfr.drop_duplicates(subset='resnum2')
                resnum1 = '-'.join([str(item) for item in dfr2.resnum2.tolist()])
                resnum2 = '-'.join([str(item) for item in dfr1.resnum1.tolist()])
                resnum1s.append(resnum1)
                resnum2s.append(resnum2)
                counter += 1
                outdf = df[df.pdbchainpair1 == inpdbid]
                # print(outdf)
                outdf['resnum1'] = resnum1
                outdf['resnum2'] = resnum2
                outdfs.append(outdf)
                print(counter)
    outdf2 = pd.concat(outdfs)
    outname = infile.split('.')[0] + '_resnum.csv'
    print(outname)
    outdf2.to_csv(outname, index=False)



def gap_in_seq():
    '''
    insert gaps in within paratope/epitope sequences
    :return:
    '''
    infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv'
    df = pd.read_csv(infile)
    print(df.shape)
    # sys.exit()
    abgapmotifs = []
    abgapmotifs2 = []
    abgapmotifs3 = []
    aggapmotifs = []
    aggapmotifs2 = []
    aggapmotifs3 = []
    for i, row in df.iterrows():
        abgaps = row.gap_pattern1.split('X')
        aggaps = row.gap_pattern2.split('X')
        pchars = list(row.sequence1)
        echars = list(row.sequence2)
        print(abgaps, aggaps, pchars, echars)
        # sys.exit()
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
    print(df.head())
    outname = infile.split('.')[0] + '_phil.csv'
    print(outname)
    # sys.exit()
    df.to_csv(outname, index=False)



def resgapmotif_dataset(infile, datasetdir, para, epi):
    '''
    prepare dl dataset for residue-gap encodeing (aggregate encoding)
    drop single chars
    :return:
    '''
    # infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil_pc.csv'
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
    print(outname, outname2)
    # sys.exit()
    paras = getattr(df, para)
    epis = getattr(df, epi)
    outcontent = '\n'.join(['\t'.join(item) for item in zip(paras, epis)])
    outcontent2 = '\n'.join(['\t'.join(item) for item in zip(epis, paras)])
    outfile  = open(outname, 'w')
    outfile.write(outcontent)
    outfile2 = open(outname2, 'w')
    outfile2.write(outcontent2)

def get_ppi_interacting_segment():
    '''
    get ppi interacting segments from pdb files in
    :return:
    '''
    pdbdir = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/3did/pdbs'
    infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_resnum.csv'
    df = pd.read_csv(infile)
    pdbchainpairs = df.pdbchainpair1.tolist()[:]
    # pdbchainpairs = [item for item in pdbchainpairs if item == '1VXQ_NE']
    total_size = len(pdbchainpairs)
    chunks = np.array_split(pdbchainpairs,4)
    ## scatter with mpi
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        data = chunks
    else:
        data = None

    data = comm.scatter(data, root=0)
    total_size = len(data)
    ## end MPI
    # for i,pdbchainpair in enumerate(pdbchainpairs):
    except_counter = 0
    for i,pdbchainpair in enumerate(data):
        try:
            parts  = pdbchainpair.split('_')
            pdbid = parts[0]
            chain1 = parts[1][0]
            chain2 =  parts[1][1]
            tag = 'inter'
            if chain1 == chain2:
                tag = 'intra'
            cdf = df[df.pdbchainpair1 == pdbchainpair]
            resnum1 = [int(item) for item in cdf.resnum1.iloc[0].split('-')]
            resnum2 = [int(item) for item in cdf.resnum2.iloc[0].split('-')]
            min1, max1 =  min(resnum1)-1, max(resnum1)+1
            min2, max2 =  min(resnum2)-1, max(resnum2)+1
            pdbfile = pdbdir + '/' + pdbid + '.pdb'
            pdbchainfile1 = pdbdir + '/' + pdbid + '_' + chain1 + '.pdb'
            pdbchainfile2 = pdbdir + '/' + pdbid + '_' + chain2 + '.pdb'
            command1 = 'pdb_selchain -%s %s > %s' % (chain1, pdbfile, pdbchainfile1)
            command2 = 'pdb_selchain -%s %s > %s' % (chain2, pdbfile, pdbchainfile2)
            os.system(command1)
            os.system(command2)
            segmentdir = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/3did/interacting_residues'
            sfile1 = segmentdir + '/' + pdbid + '_' + chain1 + '_%s%s_' % (chain1, chain2) + '_interacting_segment_%s.pdb' % tag
            sfile2 = segmentdir + '/' + pdbid + '_' + chain2 + '_%s%s_' % (chain1, chain2) + '_interacting_segment_%s.pdb' % tag
            scommand1 = 'pdb_selres -%s:%s %s > %s' % (min1, max1, pdbchainfile1, sfile1)
            scommand2 = 'pdb_selres -%s:%s %s > %s' % (min2, max2, pdbchainfile2, sfile2)
            os.system(scommand1)
            os.system(scommand2)
            ## remove header and other bits
            content1 = open(sfile1).read().splitlines()
            content1 = '\n'.join([item for item in content1 if item.startswith('ATOM')])
            outfile1 = open(sfile1, 'w')
            outfile1.write(content1)
            content2 = open(sfile2).read().splitlines()
            content2 = '\n'.join([item for item in content2 if item.startswith('ATOM')])
            outfile2 = open(sfile2, 'w')
            outfile2.write(content2)
            ## end remove header and other bits
            print('str %s done of %s, rank %s' %  (i, total_size, rank))
            # print('str %s of %s' % (i, total_size))
            # print(cdf)
            # sys.exit()
        except:
            except_counter += 1
            print(except_counter)
            print('excepting...')

def check_structures():
    '''
    check 0 byte structure
    :return:
    '''
    egstr = '5ZR1_AA'
    infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_resnum.csv'
    df = pd.read_csv(infile)
    sdf = df[df.pdbchainpair2 ==egstr]
    print(sdf)

#run stuff
# sl_dl_summary()
# get_ppi_resnum()
# gap_in_seq()
# infile = 'abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_phil.csv'
# resgapmotif_dataset(infile, 'ppiressingle', 'abgapmotif3', 'aggapmotif3')
# lost the codes for downloading pdbs (with mpi) due to rm error.
get_ppi_interacting_segment()
# check_structures()



