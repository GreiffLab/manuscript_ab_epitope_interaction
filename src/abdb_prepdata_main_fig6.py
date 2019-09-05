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


def randomize_input_files():
    '''
    randomize input files for the dl bit
    :return:
    '''
    infiles = fifi('../dl/dataset', '.tsv')
    print(len(infiles))
    infiles = [item for item in infiles if 'motif_epipara.tsv' not in item]
    for infile in infiles:
        contents = open(infile).readlines()
        print(contents[:5])
        sources = [item.split('\t')[0] for item in contents]
        targets = [item.split('\t')[1] for item in contents]
        random.shuffle(sources)
        outpairs  = ['\t'.join([source, target]) + '\n' for source, target in zip(sources, targets)]
        print(outpairs[:5])
        print(contents[:5])

        sys.exit()


# run stuff
# make_tabsep_paraepi()
# make_motif_epipara_file()
# abstract_gap()
# char_position('../dl/dataset/motif_epiparadash.tsv')
# make_paraepi_file()
# char_position('../dl/dataset/motif_paraepidash.tsv')
# make_seq_epipara_file()
randomize_input_files()




