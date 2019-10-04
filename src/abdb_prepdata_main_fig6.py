# note to self: this fig is the dl bit, files are dumped in a different directory (../dl)
# output all preped datasets to ../dl/dataset
# train model on FRAM
# remember gaps are abstrated further (to simply dash - )

# preparations and ajdusments to the dataset

import pandas as pd
import sys
import os
from find_files import find_files as fifi
import numpy as np
import random
import jellyfish

def make_tabsep_paraepi():
    '''
    make tab separated para epi file
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber.csv'
    df = pd.read_csv(infile)
    df = df.dropna(subset=['paratope', 'epitope', 'ab_motif', 'ag_motif'])
    print(df.head())
    outcontent = ''
    for i, row in df.iterrows():
        # print(row.epitope, row.paratope)
        content = '\t'.join([row.paratope, row.epitope]) + '\n'
        outcontent += content
    outname = '../dl/dataset/paraepi.tsv'
    outfile = open(outname, 'w')
    outfile.write(outcontent)
    print(outcontent)

def make_motif_epipara_file():
    '''
    make motif infile for the deep learning bit
    :return:
    '''
    # infile = '/Users/rahmadakbar/greifflab/aims/aimugen/src/abdb_outfiles/respairs_segment_notationx_merged_angle_len.csv'
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber.csv'
    df = pd.read_csv(infile)
    print(df.shape)
    print(df.columns)
    # df = df.dropna(subset=['ab_motif', 'ag_motif'])
    df = df.dropna(subset=['paratope', 'epitope', 'ab_motif', 'ag_motif'])
    outcontent = ''
    for i, row in df.iterrows():
        agmotif = row.ag_motif
        abmotif = row.ab_motif
        # print(agmotif, abmotif)
        abparts = abmotif.split('X')
        agparts  = agmotif.split('X')
        abcontents = []
        for abpart in abparts[:-1]:
            if abpart == '':
                abcontents.append('X')
            else:
                abcontents.append(abpart)
                abcontents.append('X')
        agcontents = []
        for agpart in agparts[:-1]:
            if agpart == '':
                agcontents.append('X')
            else:
                agcontents.append(agpart)
                agcontents.append('X')

        # print(abcontents)
        # print(agcontents)
        content = '\t'.join([' '.join(agcontents), ' '.join(abcontents)]) + '\n'
        # print([content])
        outcontent += content
    outpath = '../dl/dataset/motif_epipara.tsv'
    outfile = open(outpath, 'w')
    outfile.write(outcontent)



def make_motif_epipara_content(input, output):
    '''
    make motif infile for the deep learning bit
    :return:
    '''
    # print(agmotif, abmotif)
    abparts = input.split('X')
    agparts  = output.split('X')
    abcontents = []
    for abpart in abparts[:-1]:
        if abpart == '':
            abcontents.append('X')
        else:
            abcontents.append(abpart)
            abcontents.append('X')
    agcontents = []
    for agpart in agparts[:-1]:
        if agpart == '':
            agcontents.append('X')
        else:
            agcontents.append(agpart)
            agcontents.append('X')

    # print(abcontents)
    # print(agcontents)
    content = '\t'.join([' '.join(agcontents), ' '.join(abcontents)])
    return content

def abstract_gap():
    '''
    abtract gaps into -
    :return:
    '''
    infile = '../dl/dataset/motif_epipara.tsv'
    contents = open(infile).read().splitlines()
    outcontent = ''
    for content in contents[:]:
        parts = content.split('\t')
        datum  = []
        for part in parts:
            chars = part.split(' ')
            chargap = []
            for char in chars:
                if char == 'X':
                    chargap.append(char)
                else:
                    chargap.append('-')
            charstr = ' '.join(chargap)
            datum.append(charstr)
        outcontent += '\t'.join(datum) + '\n'
    outname = '../dl/dataset/motif_epiparadash.tsv'
    outfile = open(outname, 'w')
    outfile.write(outcontent)

def abstract_gap_content(content):
    '''
    abstract the gap in content
    :param content:
    :return:
    '''
    input, output = content.split('\t')
    print([content], 'hey')
    content = make_motif_epipara_content(input, output)
    print([content])
    parts = content.split('\t')
    datum  = []
    for part in parts:
        chars = part.split(' ')
        chargap = []
        for char in chars:
            if char == 'X':
                chargap.append(char)
            else:
                chargap.append('-')
        charstr = ' '.join(chargap)
        datum.append(charstr)
    outcontent = '\t'.join(datum) + '\n'
    return outcontent

def char_position(infile):
    '''
    adds position to the character. eg. XXX > X1 X2 X3, X-X, X1 -2 X3
    :return:
    '''
    # infile = '../dl/dataset/motif_epiparadash.tsv'
    contents = open(infile).read().splitlines()
    outcontent = ''
    for content in contents[:]:
        parts = content.split('\t')
        datum  = []
        for part in parts:
            part = [item + str(i+1) for i,item in enumerate(part.split()) ]
            part = ' '.join(part)
            datum.append(part)
        outline = '\t'.join(datum) + '\n'
        outcontent += outline
    # outname = '../dl/dataset/motif_epiparadash_pos.tsv'
    outname = infile[:-4] + '_pos.tsv'
    print(outname)
    outfile = open(outname , 'w')
    outfile.write(outcontent)


def char_content(content):
    '''
    adds position to the character. eg. XXX > X1 X2 X3, X-X, X1 -2 X3
    :return:
    '''
    parts = content.split('\t')
    datum  = []
    for part in parts:
        part = [item + str(i+1) for i,item in enumerate(part.split()) ]
        part = ' '.join(part)
        datum.append(part)
    outline = '\t'.join(datum) + '\n'
    return outline


def make_paraepi_file():
    '''
    flip the input file, now train from para to epi
    :return:
    '''
    infile = '../dl/dataset/motif_epiparadash.tsv'
    contents = open(infile).read().splitlines()
    outcontent = ''
    for content in contents[:]:
        epi, para = content.split('\t')
        newcontent = '\t'.join([para, epi]) + '\n'
        outcontent += newcontent
    infile_parts = infile.split('_')
    if 'pos' in infile:
        outname = '../dl/dataset/motif_paraepidash_pos.tsv'
    else:
        outname = '../dl/dataset/motif_paraepidash.tsv'
    outfile = open(outname, 'w')
    outfile.write(outcontent)


def make_seq_epipara_file():
    '''
    flip the input file, now train from para to epi
    :return:
    '''
    infile = '../dl/dataset/paraepi.tsv'
    contents = open(infile).read().splitlines()
    outcontent = ''
    for content in contents[:]:
        para, epi = content.split('\t')
        newcontent = '\t'.join([epi, para]) + '\n'
        outcontent += newcontent
    outname = '../dl/dataset/epipara.tsv'
    outfile = open(outname, 'w')
    outfile.write(outcontent)


def make_tabsep_ppi(infile, dirtag, seqinput, seqoutput, motifinput, motifoutput):
    '''
    make tab separated para epi file
    :return:
    '''
    df = pd.read_csv(infile)
    print(df)
    print(df.info())
    print(df.shape)
    df = df.dropna(subset=[seqinput, seqoutput, motifinput, motifoutput])
    print(df.shape)
    print(df.head())
    outcontent1 = ''
    outcontent2 = ''
    outcontent3 = ''
    outcontent4 = ''
    outcontent5 = ''
    outcontent6 = ''
    for i, row in df.iterrows():
        # print(row.epitope, row.paratope)
        content1 = '\t'.join([getattr(row, seqinput), getattr(row, seqoutput)]) + '\n'
        outcontent1 += content1
        content2 = '\t'.join([getattr(row, seqoutput), getattr(row, seqinput)]) + '\n'
        outcontent2 += content2
        content3 = '\t'.join([getattr(row, motifinput), getattr(row, motifoutput)])
        content3 = abstract_gap_content(content3)
        outcontent3 += content3
        content4 = '\t'.join([getattr(row, motifoutput), getattr(row, motifinput)])
        content4 = abstract_gap_content(content4)
        outcontent4 += content4
        content5 = char_content(content3)
        outcontent5 += content5
        content6 = char_content(content4)
        outcontent6 += content6
    print(outcontent6)
    items  = [item.split('\t') for item in outcontent6.splitlines()]
    print(items)
    df = pd.DataFrame(items, columns=['m1', 'm2'])
    print(df.head())
    print(df.m1.unique().shape)
    print(df.m2.unique().shape)
    print(df.shape)
    outcontents = [outcontent1, outcontent2, outcontent3, outcontent4, outcontent5, outcontent6]
    outtags  = ['paraepi', 'epipara', 'motif_paraepidash', 'motif_epiparadash', 'motif_paraepidash_pos',
                'motif_epiparadash_pos']
    for outcontent, outtag in zip(outcontents, outtags):
        outname = '../dl/%s/%s.tsv' % (dirtag,outtag)
        print(outname)
        outfile = open(outname, 'w')
        outfile.write(outcontent)




# run stuff
# make_tabsep_paraepi()
# make_motif_epipara_file()
# abstract_gap()
# char_position('../dl/dataset/motif_epiparadash.tsv')
# make_paraepi_file()
# char_position('../dl/dataset/motif_paraepidash.tsv')
# make_seq_epipara_file()

### new prep files starts here
make_tabsep_ppi('abdb_outfiles_2019/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired.csv',
            'dataset_ppi','sequence1', 'sequence2', 'gap_pattern1', 'gap_pattern2')

# make_tabsep_ppi('abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber.csv',
#                 'dataset','paratope', 'epitope', 'ab_motif', 'ag_motif')


