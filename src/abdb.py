# import stuff
import find_files as fifi
import numpy as np
import sys
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
import os
import pandas as pd
from pandas.api.types import CategoricalDtype
from Bio import PDB
from Bio.Seq import Seq
import string
import natsort
import random
from plotnine import *
from plotnine.data import *
import re
import networkx as nx
import jellyfish
import math
import itertools
import scipy


# set general style for plots
sns.set(style='ticks', palette = 'pastel')

#sets df to display all columns
pd.set_option('display.max_column', None)

def view_temp():
    '''
    View temporary plot
    '''
    sns.despine()
    temppdf = 'abdb_figures/temp.pdf'
    plt.savefig(temppdf)
    os.system('open %s' % temppdf)

def gview_temp(g):
    '''
    View temporary plot. uses ggplot via plotnine
    '''
    temppdf = 'abdb_figures/temp.pdf'
    g = (g
         # + theme_classic()
         + theme(strip_background=element_rect(fill='#2266FF', size=1.4, alpha=.30),
            panel_background=element_rect(fill='#2266FF', alpha=.05))
         )
    g.save(temppdf)
    os.system('open %s' % temppdf)

def save_pdf(filename):
    '''
    View temporary plot
    '''
    sns.despine()
    outfile = 'abdb_figures/%s.pdf' % filename
    plt.savefig(outfile)
    os.system('open %s' % outfile)


def gsave_pdf(g, filename):
    '''
    View temporary plot. uses ggplot via plotnine
    '''
    outfile = 'abdb_figures/%s.pdf' % filename
    g = (g + theme(strip_background=element_rect(fill='#2266FF', size=1.4, alpha=.30),
                   panel_background=element_rect(fill='#2266FF', alpha=.05))
         )
    g.save(outfile)
    os.system('open %s' % outfile)


def save_png(filename):
    '''
    View temporary plot
    '''
    sns.despine()
    outfile = 'abdb_figures/%s.png' % filename
    plt.savefig(outfile)
    os.system('open %s' % outfile)


def get_single_antigens():
    '''
    Get files with single antigen
    '''
    file_paths = fifi.find_files('../datasets/NR_LH_Protein_Martin','.pdb')
    seqres_index = []
    singles = []
    for fp in file_paths:
        chains = []
        contents = open(fp).readlines()
        newcontents = ''
        hetflag = False
        for i, line in enumerate(contents):
            if line.startswith('REMARK 950 CHAIN '):
                chain = line[11:18]
                chain_char = chain[-1]
                chain_orig_char = line[-1]
                # if chain_char != chain_orig_char:
                #     print(fp, chain_char, chain_orig_char)
                chains.append(chain)
            if line.startswith('HETATM') == False:
                newcontents += line
            elif line.startswith('HETATM') == True:
                hetflag =  True
        nchains = len(set(chains))
        if nchains == 3:
            outfile = open(fp, 'w')
            outfile.write(newcontents)
            if hetflag:
                print('file %s has 3 chains and hetatm, removing hetatm, adding it to the dataset' % fp)
            else:
                print('file %s has 3 chains, adding it to the dataset' % fp)
            singles.append(fp)
    return singles


def get_residue_pairs():
    '''
    Gets antibody-antigen residue pairs according to some distance d.
    Considers only CA atoms.
    '''
    file_paths = get_single_antigens()
    top15s = 'distance,agresname,agchain,agresnum,abresname,abchain,abresnum,name\n'
    for fp in file_paths[:]:
        parts = open(fp).read().split('TER  ')
        antigen = parts[-2]
        antigen = [line for line in antigen.split('\n') if 'CA' in line[12:16] and line.startswith('ATOM')]
        antibody = parts[0] + parts[1]
        antibody = [line for line in antibody.split('\n') if 'CA' in line[12:16] and line.startswith('ATOM')]
        min_pairs = []
        sample_name = fp.split('/')[-1].split('.')[0]
        for ag in antigen:
            min_pair = [] # nearest pair
            ag_c = np.array([ag[30:38], ag[38:46], ag[46:54]], dtype='float16')
            ag_inf = [ag[17:20], ag[21], ag[22:26]]
            min_d = 1000
            for ab in antibody:
                ab_c = np.array([ab[30:38], ab[38:46], ab[46:54]], dtype='float16')
                d = np.sqrt(np.sum((ag_c-ab_c)**2))
                if d < min_d:
                    min_d = [d]
                    ab_inf = [ab[17:20], ab[21], ab[22:26]]
                    min_pair = [d] + ag_inf + ab_inf + [sample_name]
            min_pairs.append(min_pair)
        top15 = sorted(sorted(min_pairs)[:15], key = lambda item: item[3])
        joined = ''
        for t in top15:
            t = ','.join([str(item).strip() for item in t]) + '\n'
            joined += t
        #joined += 'TER\n'
        top15s += joined
    outpath = 'abdb_outfiles/respair.csv'
    outfile = open(outpath, 'w')
    outfile.write(top15s)
    outfile.close()


def get_residue_pairs_ab():
    '''
    Gets antibody-antigen residue pairs according to some distance d.
    Considers only CA atoms.
    iterates over antibodies
    '''
    file_paths = get_single_antigens()
    top15s = 'distance,abresname,abchain,abresnum,agresname,agchain,agresnum,name\n'
    for fp in file_paths[:]:
        parts = open(fp).read().split('TER  ')
        antigen = parts[-2]
        antigen = [line for line in antigen.split('\n') if 'CA' in line[12:16] and line.startswith('ATOM')]
        antibody = parts[0] + parts[1]
        antibody = [line for line in antibody.split('\n') if 'CA' in line[12:16] and line.startswith('ATOM')]
        min_pairs = []
        sample_name = fp.split('/')[-1].split('.')[0]
        for ag in antibody: #use antibody here
            min_pair = [] # nearest pair
            ag_c = np.array([ag[30:38], ag[38:46], ag[46:54]], dtype='float16')
            ag_inf = [ag[17:20], ag[21], ag[22:26]]
            min_d = 1000
            for ab in antigen: # use antigen here
                ab_c = np.array([ab[30:38], ab[38:46], ab[46:54]], dtype='float16')
                d = np.sqrt(np.sum((ag_c-ab_c)**2))
                if d < min_d:
                    min_d = [d]
                    ab_inf = [ab[17:20], ab[21], ab[22:26]]
                    min_pair = [d] + ag_inf + ab_inf + [sample_name]
            min_pairs.append(min_pair)
        top15 = sorted(sorted(min_pairs)[:15], key = lambda item: item[3])
        joined = ''
        for t in top15:
            t = ','.join([str(item).strip() for item in t]) + '\n'
            joined += t
        #joined += 'TER\n'
        top15s += joined
    outpath = 'abdb_outfiles/respair_ab_top15.csv'
    outfile = open(outpath, 'w')
    outfile.write(top15s)
    outfile.close()


def get_residue_pairs_ab2(file_paths, outpath, cutoff):
    '''
    Gets antibody-antigen residue pairs according to some distance d.
    Considers any atoms.
    sorts antibodies and antigens
    Uses Biopython because I'm lazy
    '''
    #file_paths = get_single_antigens()
    chain_lens = []
    pairs_absort = []
    pairs_agsort = []
    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    aalist= aadict.keys()
    for fp in file_paths[:]:
        cutoff = cutoff
        parser = PDB.PDBParser(PERMISSIVE=1)
        structure_name = fp.split('/')[-1].split('.')[0]
        structure = parser.get_structure(structure_name, fp)
        chains = structure.get_chains()
        chain_list = PDB.Selection.unfold_entities(structure,'C')
        chain_lens.append(len(chain_list))
        print(chain_list)
        abchains = ['L','H']
        for i, abchain in enumerate(abchains):
            #chain order follows the pdb file. antigen always last,hence, chain_list[-1]
            atoms = PDB.Selection.unfold_entities(chain_list[i], 'A') + PDB.Selection.unfold_entities(chain_list[-1],'A')
            searcher = PDB.NeighborSearch(atoms)
            atom_pairs = searcher.search_all(radius=cutoff,level='A')
            hetero_pairs = []
            for ap in atom_pairs:
                abres, agres = ap[0].get_parent(), ap[1].get_parent()
                abchain, agchain = abres.get_parent(), agres.get_parent()
                #print('tag', chain1, chain2)
                if abchain != agchain:
                    hetero_pairs.append((abres, agres))
                    #print(dir(abres), abchain, agres.get_full_id(), agchain)
                    #break
            print(len(hetero_pairs), len(set(hetero_pairs)))
            respairs = []
            for hp in set(hetero_pairs):
                #print(help(hp[0].get_full_id()))
                name = hp[0].get_full_id()[0]
                abres,agres = hp[0].resname, hp[1].resname
                abresnum, agresnum = hp[0].get_full_id()[3][1], hp[1].get_full_id()[3][1]
                #abresnumi = abres + ''.join([str(item) for item in hp[0].get_full_id()[3]]).strip()
                abresnumi = ''.join([str(item) for item in hp[0].get_full_id()[3]]).strip()
                #agresnumi = agres + ''.join([str(item) for item in hp[1].get_full_id()[3]]).strip()
                agresnumi = ''.join([str(item) for item in hp[1].get_full_id()[3]]).strip()
                abchain, agchain = hp[0].get_full_id()[2], hp[1].get_full_id()[2]
                respair = [name, abres, abresnum, abresnumi, abchain, agres, agresnum, agresnumi, agchain]
                if respair not in respairs: respairs.append(respair)
            #if name == '4LSS_1' and abchain=='H':
                #print(pd.DataFrame(respairs))
                #print(pd.DataFrame(natsort.natsorted(respairs, key= lambda item: item[3])))
                #sys.exit()
            pairs_absort = pairs_absort + natsort.natsorted(respairs, key= lambda item: item[3])
            pairs_agsort= pairs_agsort + sorted(respairs, key= lambda item: item[-2])
    columns = ['pdbid', 'abres', 'abresnum', 'abresnumi', 'abchain', 'agres', 'agresnum', 'agresnumi', 'agchain']
    df = pd.DataFrame(pairs_absort, columns = columns)
    df2 = pd.DataFrame(pairs_agsort, columns = columns)
    aboutname = outpath + '/respairs_absort_cutoff%s.csv' % cutoff
    agoutname = outpath + '/respairs_agsort_cutoff%s.csv' % cutoff
    df2.to_csv(agoutname, index=False)
    df.to_csv(aboutname, index=False)
    return aboutname, agoutname

def add_gaps():
    '''
    adds gaps
    '''
    df, df2 = pd.read_csv('abdb_outfiles/respairs_absort.csv'), pd.read_csv('abdb_outfiles/respairs_agsort.csv')
    abshifts = []
    agshifts = []
    for i in range(df.shape[0]-1):
        sample,sample2 = df.iloc[i], df2.iloc[i]
        next_sample, next_sample2 = df.iloc[i+1], df2.iloc[i+1]
        if sample.pdbid == next_sample.pdbid and sample.abchain == next_sample.abchain:
            abshift = next_sample.abresnum - sample.abresnum
            abshifts.append(abshift)
            agshift = next_sample2.agresnum - sample2.agresnum
            agshifts.append(agshift)
        else:
            abshifts.append(0)
            agshifts.append(0)
    print(df.shape)
    print(len(abshifts))
    df['abshift'] = abshifts + [0]
    df['agshift'] = agshifts + [0]
    df.to_csv('abdb_outfiles/respairs_gaps.csv',index=False)

def get_unique_abresnumi(infile, outpath):
    '''
    get unique interacting residues. account for inserted residues
    '''
    #infile =  'abdb_outfiles/respairs_absort.csv'
    df = pd.read_csv(infile)#, pd.read_csv('abdb_outfiles/respairs_agsort.csv')
    dfs = []
    chains = ['L', 'H']
    for pdbid in df.pdbid.unique():
        pdbdf = df[df.pdbid == pdbid]
        for chain in chains:
            chaindf = pdbdf[pdbdf.abchain == chain].drop_duplicates(subset='abresnumi')
            dfs.append(chaindf)
    newdf = pd.concat(dfs)
    outname = outpath + '/' + infile.split('/')[-1].split('.')[0] + '_abresnumi.csv'
    print(outname)
    print(newdf.shape, df.shape)
    newdf.to_csv(outname, index=False)
    return outname

def add_segments(infile, outpath):
    '''
    adds CDR and FR segments.
    '''
    #infile = 'abdb_outfiles/respairs_absort_abresnumi.csv'
    df = pd.read_csv(infile)
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
    chains = ['L', 'H']
    segment_dicts = dict([('L',lsegment_dict), ('H',hsegment_dict)])
    segment_labels = []
    segment_labels_frcdr = []
    print(df.head())
    for i, row in df.iterrows():
        segment_dict = segment_dicts[row.abchain]
        segment_label = segment_dict[row.abresnum]
        segment_labels.append(segment_label)
        frcdr = '%s%s'% (row.abchain,segment_label[-1])
        segment_labels_frcdr.append(frcdr)
    df['loop'] = segment_labels_frcdr
    df['segment'] = segment_labels
    outname = outpath + '/' + infile.split('/')[-1].split('.')[0] + '_segments.csv'
    df.to_csv(outname, index=False)
    return outname


def add_abshift(infile,outpath):
    '''
    Add positional shifts based on residue numbers, account for inserted residues and segments.
    eg line i position 30, line i+1 = 32, shifted position = 2
    eg line i position 30, line i+1 = 32A, shifted position = 3 (Martin insertion position)
    '''
    #infile = 'abdb_outfiles/respairs_absort_abresnumi_segments.csv'
    df = pd.read_csv(infile)
    abshifts = []
    alphabets = list(string.ascii_uppercase)
    alphabets = dict([item,i] for i,item in enumerate(alphabets))
    for i in range(len(df)-1)[:]:
        current_row = df.iloc[i]
        next_row =  df.iloc[i+1]
        if current_row.segment == next_row.segment and current_row.pdbid == next_row.pdbid:
            try:
                abshift = int(next_row.abresnumi) - int(current_row.abresnumi)
            except:
                print('Found inserted residues, need to do some gymnastic...')
                if current_row.abresnum != next_row.abresnum:
                    abshift = next_row.abresnum - current_row.abresnum
                else:
                    print('more gymnastic...')
                    if len(current_row.abresnumi) != len(next_row.abresnumi):
                        abshift = alphabets[next_row.abresnumi[-1]] + 1
                    else:
                        abshift = alphabets[next_row.abresnumi[-1]] - alphabets[current_row.abresnumi[-1]]
        else:
            abshift = np.nan
        abshifts.append(abshift)
    abshifts.append(np.nan) # end piece
    df['abshift'] = abshifts
    outname = outpath + '/' + infile.split('/')[-1].split('.')[0] + '_abshift.csv'
    df.to_csv(outname, index=False)
    return outname


def add_abshiftl(infile,outpath):
    '''
    Add shifts based on residue numbers, account for inserted residues and loop.
    '''
    #infile = 'abdb_outfiles/respairs_absort_abresnumi_segments_abshift.csv'
    df = pd.read_csv(infile)
    abshifts = []
    alphabets = list(string.ascii_uppercase)
    alphabets = dict([item,i] for i,item in enumerate(alphabets))
    for i in range(len(df)-1)[:]:
        current_row = df.iloc[i]
        next_row =  df.iloc[i+1]
        if current_row.loop == next_row.loop and current_row.pdbid == next_row.pdbid:
            try:
                abshift = int(next_row.abresnumi) - int(current_row.abresnumi)
            except:
                print('Found inserted residues, need to do some gymnastic...')
                if current_row.abresnum != next_row.abresnum:
                    abshift = next_row.abresnum - current_row.abresnum
                else:
                    print('more gymnastic...')
                    if len(current_row.abresnumi) != len(next_row.abresnumi):
                        abshift = alphabets[next_row.abresnumi[-1]] + 1
                    else:
                        abshift = alphabets[next_row.abresnumi[-1]] - alphabets[current_row.abresnumi[-1]]
        else:
            abshift = np.nan
        abshifts.append(abshift)
    abshifts.append(np.nan) # end piece
    df['abshiftl'] = abshifts
    print(df.head(50))
    outname = outpath + '/' + infile.split('/')[-1].split('.')[0] + '_abshiftl.csv'
    df.to_csv(outname, index=False)
    return outname


def lh_counts():
    '''
    Checks light and heavy chain counts
    '''
    samples = open('abdb_outfiles/respair.csv').read().split('TER')
    H = 0
    L = 0
    for i, sample in enumerate(samples):
        residues = sample.split()
        for r in residues:
            r = r.split(',')
            if r[5] == 'H':
                H += 1
            elif r[5] == 'L':
                L += 1
    hfrac = float(H)/sum([H,L])
    lfrac = float(L)/sum([H,L])
    x = ['H', 'L']
    y = [H,L]
    sns.barplot(x,y)
    plt.title('Residues in heavy and light chains')
    plt.xlabel('Chain')
    plt.ylabel('Count')
    filename = 'abdb_figures/hl.pdf'
    plt.savefig(filename)
    os.system('open %s' % filename)
    #print(help(sns.barplot))
    return H, L


def respair_exams(filepath):
    '''
    Examines respair.csv
    obs: 1/3 of epitope residues are in the light chain.
           spatial distance is between 3-10 angstrom.
    '''
    filename = filepath.split('/')[-1].split('.')[0]
    dataset = pd.read_csv(filepath)
    sns.set(style='ticks', palette='BuGn')
    sns.catplot(x='abchain', kind='count', data=dataset)
    outname = 'abdb_figures/%s_abchain_count.pdf' % filename
    plt.savefig(outname)
    os.system('open %s' % outname)
    outname2 = 'abdb_figures/%s_abchain_distance.pdf' % filename
    sns.catplot(x='abchain',y='distance', kind='box', data=dataset)
    plt.savefig(outname2)
    os.system('open %s' % outname2)


def positional_distance(filepath):
    '''
    Plots distance between positions
    '''
    filename = filepath.split('/')[-1].split('.')[0]
    #infile = 'abdb_outfiles/respair.csv'
    df = pd.read_csv(filepath)
    name_distance = df.xs(['name','agresnum','abresnum'],axis=1).values
    resnum_diffs = []
    for i in range(len(name_distance)-1):
        if name_distance[i,0] == name_distance[i+1,0]:
            diff = np.abs(name_distance[i,1:] - name_distance[i+1,1:])
            if diff[0] < 15:
                resnum_diffs.append(diff)
    resnum_diffs = pd.DataFrame(resnum_diffs, columns=['agdist','abdist'])
    #sns.set(style='ticks', palette = 'BuGn_d')
    sns.jointplot(x='agdist', y='abdist', data=resnum_diffs, kind='kde')
    #outname = 'abdb_figures/%s_agdist_abdist.pdf' % filename
    #plt.savefig(outname)
    #os.system('open %s' % outname)
    #print(help(df))


def get_h_epitope():
    '''
    Filters respair.csv for heavy chain epitopes.
    '''
    df = pd.read_csv('abdb_outfiles/respair.csv')
    #print(dir(df))
    filepath = 'abdb_outfiles/respair_h.csv'
    newvalues=  []
    df.to_csv(filepath)
    for line in df.values:
        if line[5]  == 'H':
            newvalues.append(line)
    newdf = pd.DataFrame(newvalues, columns=df.columns.values)
    newdf.to_csv(filepath, index=False)
    #plt.show()


def get_h_epitope_ab():
    '''
    Filters respair.csv for heavy chain epitopes.
    '''
    df = pd.read_csv('abdb_outfiles/respair_ab_full.csv')
    #print(dir(df))
    filepath = 'abdb_outfiles/respair_ab_full_h.csv'
    newvalues=  []
    df.to_csv(filepath)
    for line in df.values:
        if line[2]  == 'H':
            newvalues.append(line)
    newdf = pd.DataFrame(newvalues, columns=df.columns.values)
    newdf.to_csv(filepath, index=False)
    #plt.show()



def get_seq_pair():
    '''
    Gets antibody-antigen sequence pairs
    '''
    infile = 'abdb_outfiles/respair_ab_full_h.csv'
    df = pd.read_csv(infile)
    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    seqs = []
    samples = {}
    values = df.values
    for val in values:
        name = val[-1]
        agres = val[1]
        abres = val[4]
        print(val)
        if name not in samples:
            samples[name] = [aadict[agres], aadict[abres]]
        elif name in samples:
            samples[name][0] += aadict[agres]
            samples[name][1] += aadict[abres]
    #print(samples.items())
    todf = [[key] + val for key, val in samples.items()]
    newdf = pd.DataFrame(todf, columns=['name','abseq','agseq'])
    outfile = infile.split('.')[0] + '_seqpair.' + infile.split('.')[-1]
    newdf.to_csv(outfile, index=False)
    #print(help(pd.DataFrame))


def find_car():
    '''
    Finds CAR segment on a sequence
    '''
    infile = 'abdb_outfiles/respair_ab_full_h_seqpair.csv'
    df = pd.read_csv(infile)
    rnums = [] # CAR start position
    len_cars = [] # length of seq
    variants = [] # CAR variants
    i_cars = []
    for i,row in df.iterrows():
        n = len(row.abseq)
        for rnum in range(n-2):
            triplet = row.abseq[rnum:rnum+3]
            if triplet == 'CAR' and i not in i_cars:
                rnums.append(rnum)
                len_cars.append(n)
                i_cars.append(i)
    for i, row in df.iterrows():
        if i not in i_cars:
            n = len(row.abseq)
            for rnum in range(n-2):
                triplet = row.abseq[rnum:rnum+3]
                if triplet.startswith('CA'):
                    rnums.append(rnum)
                    len_cars.append(n)
                    i_cars.append(i)
                    variants.append(triplet)
    pltdf = pd.DataFrame({'car_start':rnums, 'seq_len': len_cars})
    print(set(variants))
    print(len(len_cars), len(rnums))
    #sns.distplot(rnums)
    sns.jointplot(x=pltdf.car_start, y=pltdf.seq_len, kind='kde')
    save_png('car_position')
    #view_temp()
    #print(help(sns.distplot))


def neighbor_gap():
    '''
    Finds the gap between residues
    '''
    infile = 'abdb_outfiles/respair_ab_full_h.csv'
    df = pd.read_csv(infile)
    #abdf = df.xs(['distance', 'abresnum', 'name'], axis=1)
    #abdf_1 = abdf.shift(-1)
    #distance_diff = abdf.distance - abdf_1.distance
    #print(abdf.distance)
    #print(help(df.xs))
    diff_df = []
    for i in range(len(df)-1)[:]:
        current_res, next_res = df[i:i+1], df[i+1:i+2]
        if current_res.name.values == next_res.name.values and current_res.abresnum.values > 80:
            cvalues, nvalues= list(current_res.values[0]), next_res.values[0]
            cvalues.append( next_res.abresnum.values[0] - current_res.abresnum.values[0])
            cvalues.append(next_res.agresnum.values[0] - current_res.agresnum.values[0])
            #if cvalues[0] < 15:
            diff_df.append(cvalues)
    headers=list(df.keys().values)
    headers.append('abdist')
    headers.append('agdist')
    diff_df = pd.DataFrame(diff_df, columns=headers)
    #sns.set(style='ticks', palette = 'BuGn_d')
    sns.jointplot(x='agdist', y='abdist', data=diff_df)#, kind='kde')
    view_temp()
    outname = 'abdb_outfiles/respair_abagdist.csv'
    diff_df.to_csv(outname, index=False)


def get_abagdist():
    '''
    Adds gap columns (abdist, agdist) to infile
    '''
    infile = 'abdb_outfiles/respair_ab_full.csv'
    df = pd.read_csv(infile)
    dfh = df[df['abchain'] == 'H']
    dfl = df[df['abchain'] == 'L']
    tags, dfs = ['H', 'L'], [dfh,dfl]
    for tag, df in zip(tags,dfs):
        diff_df = []
        for i in range(len(df)-1)[:]:
            current_res, next_res = df[i:i+1], df[i+1:i+2]
            if current_res.name.values == next_res.name.values and current_res.abresnum.values > 0:
                cvalues, nvalues= list(current_res.values[0]), next_res.values[0]
                cvalues.append(next_res.abresnum.values[0] - current_res.abresnum.values[0])
                cvalues.append(next_res.agresnum.values[0] - current_res.agresnum.values[0])
                if cvalues[0] < 15:
                    diff_df.append(cvalues)
        headers=list(df.keys().values)
        headers.append('abdist')
        headers.append('agdist')
        diff_df = pd.DataFrame(diff_df, columns=headers)
        print(diff_df.head(20))
        #sns.jointplot(x='agdist', y='abdist', data=diff_df)#, kind='kde')
        #view_temp()
        outname = infile.split('.')[0] + '_%s_abagdist.csv' % tag
        diff_df.to_csv(outname, index=False)


def exam_neighbor_gap():
    '''
    Examines neighbor gap
    '''
    infile = 'abdb_outfiles/respair_abagdist.csv'
    df = pd.read_csv(infile)
    df2 = df.loc[df['distance'] < 10]
    #plt.subplot(2,1,1)
    cols = ['abdist', 'agdist']
    f, axes = plt.subplots(len(cols),1)
    for i,col in enumerate(cols):
        sns.violinplot(x = col, data=df2, ax=axes[i])
    plt.tight_layout()
    #view_temp()
    save_png('neighbor_gap')
    df3 = df2.loc[np.abs(df2['agdist']) == 0]
    print(df3.head(20))
    print(df3.shape)


def lab180813_distribution():
    '''
    add residue number distributionn on top of the gap plot (antibody and antigen)
    look into the light chain as well
    distribution of interacting residues (light vs chain)
    communicate with Jeliazko
    check linearity (resolve conformational vs linear segments)
    '''
    ### distribution
    abagdists = fifi.find_files('abdb_outfiles','abagdist')
    f, axes = plt.subplots(1,len(abagdists)*2)
    axes_i = 0
    features = ['abresname', 'agresname']
    xmin, xmax = 0, 0
    ss = []
    labels = []
    for path in abagdists:
        tag = path.split('_')[-2]
        df = pd.read_csv(path)
        for i, feat in enumerate(features):
            label = feat + ' %s' % tag
            labels.append(label)
            ax = axes[axes_i]
            #s =  np.log(df[feat].value_counts()).sort_index()
            s = df[feat].value_counts().sort_index()
            ss.append(s)
            sns.barplot(y=s.index, x=s, ax=ax, palette='pastel')
            ax.set(ylabel=label, xlabel='count', xlim=[0,5000])
            ax.set_yticklabels(labels=s.index, fontsize=8)
            axes_i += 1
    outname = 'interres_distributions'
    plt.tight_layout()
    save_png(outname)
    save_pdf(outname)
    print(len(ss))
    f,axes = plt.subplots(1,len(labels)-2)
    for i in range(len(labels)-2):
        diff = ss[i+2]-ss[i]
        ax = axes[i]
        print(i, labels)
        sns.barplot(y=diff.index, x=diff,ax=ax, palette='pastel')
        ax.set(ylabel='%s-%s'%(labels[i+2], labels[i]), xlabel='count difference')
        #ax.set_yticklabels(labels=diff.index, fontsize=8)
    plt.tight_layout()
    outname = 'distribution_difference'
    save_png(outname)
    save_pdf(outname)

def lab180813_linearity():
    abagdists = fifi.find_files('abdb_outfiles','abagdist')
    f2, axes2 = plt.subplots(1,len(abagdists)*2)
    features2 = ['abdist', 'agdist']
    axes_i = 0
    for path in abagdists:
        tag = path.split('_')[-2]
        df = pd.read_csv(path)
        for i, feat2 in enumerate(features2):
            confs = {}
            label2 = feat2 + ' %s (# of samples)' % tag
            ax2 = axes2[axes_i]
            for i in range(8):
                s2 = df[df[feat2] >= i]
                confs[i] = len(set(s2['name'].values))
            confs = pd.Series(confs)
            #print(confs)
            sns.barplot(y=confs, x=confs.index, ax=ax2)
            ax2.set_xticklabels(confs.index, fontsize=8)
            ax2.set(ylabel=label2, xlabel='gap')
            axes_i += 1
    plt.tight_layout()
    outname = 'interres_gaps'
    save_png(outname)
    save_pdf(outname)


def lab180821_pre():
    '''
    Pre lab meeting stuff
    '''
    infile = 'abdb_outfiles/respair_ab_full_h_seqpair.csv'
    print(infile)
    df = pd.read_csv(infile)
    df2 = pd.DataFrame()
    df2['paralen'] = [len(set(item)) for item in df.abseq]
    df2['epilen'] = [len(set(item)) for item in df.agseq]
    df2['paratope'] = [''.join(set(item)) for item in df.abseq]
    df2['epitope'] = [''.join(set(item)) for item in df.agseq]
    k = 3
    sparatopes = []
    sepitopes = []
    for paratope, epitope in zip(df2.paratope, df2.epitope):
        i,i2 = np.random.randint(len(paratope)), np.random.randint(len(epitope))
        sparatope = paratope[:i] + 'V'*k + paratope[i+k:]
        sepitope = epitope[:i2] + 'G'*k + epitope[i2+k:]
        sparatopes.append(sparatope)
        sepitopes.append(sepitope)
    df2['sparatope'] = sparatopes
    df2['sepitope'] = sepitopes
    print(df2.head(20))



def plot_gaps():
    '''
    check conf v linear
    '''
    infile = 'abdb_outfiles/respairs_gaps.csv'
    df = pd.read_csv(infile)
    vars = ['abshift', 'agshift']
    gapdf = pd.DataFrame()
    gapcounts = []
    for var in vars:
        hcounts = []
        lcounts = []
        for i in range(8):
            temph = df[df.abchain == 'H']
            count = set(temph[temph[var] >= i].pdbid)
            hcounts.append(len(count))
            templ = df[df.abchain == 'L']
            lcount = set(templ[templ[var] >= i].pdbid)
            lcounts.append(len(lcount))
        gapcounts.append(hcounts)
        gapcounts.append(lcounts)
    print(gapcounts)
    gapcdf = pd.DataFrame(np.array(gapcounts).T, columns = ['abshift H', 'abshift L', 'agshift H', 'agshift L'])
    f, axes = plt.subplots(1,4)
    orders = ['abshift H', 'agshift H', 'abshift L', 'agshift L']
    for i, order in enumerate(orders):
        sns.barplot(x=gapcdf[order].index, y=gapcdf[order],palette='pastel', ax=axes[i])
        ylabel = order + ' (# of samples)'
        axes[i].set(ylim=[0,800], ylabel=ylabel, xlabel='gap')
    plt.tight_layout()
    plotname = infile.split('/')[-1].split('.')[0] + '_samplespergap'
    save_png(plotname)



def plot_cdr():
    '''
    plots residue, distribution, count difference, and gaps for each CDR
    Uses the range from abresnum distribution.
    cdr1: 20 - 45
    cdr2: 46 - 70
    cdr3: > 80 (80-130)
    '''
    infile = 'abdb_outfiles/respairs_gaps_cdr.csv'
    df = pd.read_csv(infile)
    #ranges = [(0,45), (46,79), (80,130)]
    ranges = [(31,35), (50,65), (95,102)] #cdrh
    aagroups = [['ASP', 'GLU', 'LYS', 'ARG','SER', 'THR', 'GLN', 'ASN', 'HIS', 'CYS', 'ALA', 'VAL','LEU', 'ILE', 'PHE', 'TRP', 'MET','PRO','GLY','TYR'],
        ['c','c', 'c', 'c','p', 'p','p','p','p','p','h','h','h','h','h','h','h','h', 'h', 'p']]
    for ir,r in enumerate(ranges):
        dfr = df[df.abresnum >= r[0]]
        dfr = dfr[dfr.abresnum <= r[1]]
        dfrl = dfr[dfr.abchain == 'L']
        dfrh  = dfr[dfr.abchain == 'H']
        dfrs = [dfrh,dfrl]
        #for dfr in dfrs:
            #sns.countplot(dfr.abres)
        f, axes = plt.subplots(1,len(dfrs)*2)
        ix = 0
        vars = ['abres','agres']
        chains = ['H','L']
        for i, chain in enumerate(chains):
            dfr = dfrs[i]
            for var in vars:
                label = var + ' ' + chain
                ax = axes[ix]
                aagdf = pd.DataFrame({'resname':aagroups[0], 'resgroup':aagroups[1]}, index=aagroups[0]).sort_index()
                vs = pd.DataFrame({'counts':dfr[var].value_counts()})
                aagdf = pd.concat([aagdf,vs], axis=1)
                sns.barplot(y='resname', x='counts', hue='resgroup', data=aagdf, palette='RdYlBu_r',ax=ax,dodge=False)
                #sns.countplot(y=sorted(dfr[var]), ax=ax, palette='pastel')
                #sns.countplot(y=sorted(df.agres), ax=axes[ix+1])
                ax.set(ylabel= label, xlim=[0,4000])
                ix+=1
        f.set_size_inches(10,5)
        stitle = 'CDR%s' % str(ir+1)
        f.suptitle(stitle)
        plt.tight_layout()
        plotname = infile.split('/')[-1].split('.')[0] + '_distribution_' + stitle
        save_pdf(plotname)
        ### gap plots
        vars = ['abshift', 'agshift']
        gapdf = pd.DataFrame()
        gapcounts = []
        for var in vars:
            hcounts = []
            lcounts = []
            for i2 in range(8):
                #temph = dfr[dfr.abchain == 'H']
                temph = dfrh
                count = set(temph[temph[var] >= i2].pdbid)
                hcounts.append(len(count))
                #templ = dfr[dfr.abchain == 'L']
                templ = dfrl
                lcount = set(templ[templ[var] >= i2].pdbid)
                lcounts.append(len(lcount))
            gapcounts.append(hcounts)
            gapcounts.append(lcounts)
        gapcdf = pd.DataFrame(np.array(gapcounts).T, columns = ['abshift H', 'abshift L', 'agshift H', 'agshift L'])
        f2, axes2 = plt.subplots(2,2)
        orders = ['abshift H', 'agshift H', 'abshift L', 'agshift L']
        orders = np.array(orders).reshape(2,2)
        for irow, order in enumerate(orders):
            for icol, col in enumerate(order):
                ax2=axes2[irow, icol]
                #sys.exit()
                sns.barplot(x=gapcdf[col].index, y=gapcdf[col],palette='pastel', ax=ax2)
                ylabel = col + ' (# of samples)'
                ax2.set(ylim=[0,800], ylabel=ylabel, xlabel='shift')
                ax2.set_ylabel(ylabel,fontsize=10)
        f2.suptitle(stitle,y=0.98, x=0.47)
        f2.tight_layout()
        plotname = infile.split('/')[-1].split('.')[0] + '_samplespergap_' + stitle
        save_pdf(plotname)

def add_cdr_column():
    '''
    Adds a column cdr to respairs_gaps.csv
    '''
    infile = 'abdb_outfiles/respairs_gaps.csv'
    df = pd.read_csv(infile)
    ranges = [(0,35), (36,65), (66,102)] #frh+cdrh
    #ranges = [(31,35), (50,65), (95,102)] #cdrh
    df['cdr'] = df.abresnum
    for i,r in enumerate(ranges):
        mask = df.cdr.between(r[0],r[1], inclusive=True)
        df.loc[mask,'cdr'] = i+1
    g1 = df['cdr'].iloc[:-1]
    g2 = df['cdr'].iloc[1:].reset_index()
    gmask = g1 != g2.cdr
    fend = pd.Series([False], name='cdr')
    gmask = gmask.append(fend, ignore_index=True)
    df.loc[gmask,'abshift'] = 0
    df = df[df.cdr <= 3]
    df.to_csv('abdb_outfiles/respairs_gaps_cdr.csv', index=False)
    df2  = pd.read_csv('abdb_outfiles/respairs_gaps_cdr.csv')
    print(df.head(50))

def get_numgaps():
    '''
    Gets numerical descriptions from the gaps
    '''
    infile = 'abdb_outfiles/respairs_gaps_cdr.csv'
    df =  pd.read_csv(infile)
    dfl = df[df.abchain == 'L']
    dfh = df[df.abchain == 'H']
    dfs = [dfh,dfl]
    numsdfs = []
    numsdfs2 = []
    for df in dfs:
        samples = df.pdbid.unique()
        nums = []
        nums2 = []
        emps = []
        for sample in samples:
            tempdf = df[df.pdbid == sample]
            for cdr in range(1,4):
                tempdf2 = tempdf[tempdf.cdr == cdr]
                #print(tempdf2.shape)
                #print(plen)
                if len(tempdf2) == 0:
                    emps.append((sample,cdr))
                else:
                    plenresnum = tempdf2.abresnum.max() - tempdf2.abresnum.min()
                    plen = len(tempdf2.abresnum.unique())
                    maxgap = tempdf2.abshift.max()
                    gaps = [item for item in tempdf2.abshift.unique() if item > 1]
                    shiftset = '-'.join([str(item) for item in tempdf2.abshift if item > 0])
                    gapset = '-'.join([str(item-1) for item in tempdf2.abshift if item > 0])
                    ngaps = len(gaps)
                    chain = tempdf2.abchain.unique()[0]
                    if ngaps != 0:
                        sgapset = '-'.join([str(item) for item in gaps])
                        nums2.append((sample, plen, plenresnum, maxgap, ngaps, sgapset, cdr,chain, shiftset,gapset))
                        for sgap in gaps:
                            nums.append((sample, plen, plenresnum, maxgap, ngaps, sgap, cdr,chain, shiftset,gapset))
                    else:
                        sgap = 0
                        nums.append((sample, plen, plenresnum, maxgap, ngaps, sgap, cdr,chain, shiftset,gapset))
                        nums2.append((sample, plen, plenresnum, maxgap, ngaps, '-', cdr,chain, shiftset,gapset))
        print(len(emps)/float(len(nums)))
        numsdf = pd.DataFrame(nums, columns = ['name','plen', 'plenresnum', 'maxgap', 'ngaps','sgap','cdr','chain', 'shiftset','gapset'])
        numsdf2 = pd.DataFrame(nums2, columns = ['name','plen', 'plenresnum', 'maxgap', 'ngaps','sgapset','cdr','chain', 'shiftset','gapset'])
        numsdfs.append(numsdf)
        numsdfs2.append(numsdf2)
    newdf = pd.concat(numsdfs,axis=0, ignore_index=True)
    newdf2 = pd.concat(numsdfs2,axis=0, ignore_index=True)
    newdf.to_csv('abdb_outfiles/respairs_abshift_numgaps.csv', index=False)
    newdf2.to_csv('abdb_outfiles/respairs_abshift_numgaps2.csv', index=False)

def get_numgaps_segment(infile,outpath):
    '''
    Gets shifts and gaps descriptions.
    '''
    #infile = 'abdb_outfiles/respairs_absort_abresnumi_segments_abshift_abshiftl.csv'
    df = pd.read_csv(infile)
    chains = ['L', 'H']
    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    segment_descs = []
    for pdbid in df.pdbid.unique():
        pdbdf = df[df.pdbid == pdbid]
        for chain in chains:
            chaindf = pdbdf[pdbdf.abchain == chain]
            segments = chaindf.segment.unique()
            for segment in segments:
                segmentdf = chaindf[chaindf.segment == segment]
                paratope = ''.join([aadict[item] for item in segmentdf.abres])
                abresnumiset = '-'.join(str(item) for item in segmentdf.abresnumi)
                plen = len(paratope)
                shiftset = '-'.join([str(int(item)) for item in segmentdf.abshift.dropna()])
                gapset = '-'.join([str(int(item)-1) for item in segmentdf.abshift.dropna()])
                if len(shiftset) == 0:
                    shiftset = str(np.nan)
                    gapset = str(np.nan)
                segment_desc = [pdbid,chain,segment,paratope, plen, shiftset, gapset, abresnumiset]
                segment_descs.append(segment_desc)
    columns = ['pdbid', 'abchain', 'segment', 'paratope', 'plen', 'shiftset', 'gapset', 'abresnumiset']
    segment_descdf = pd.DataFrame(segment_descs, columns = columns)
    inname = infile.split('/')[-1].split('.')[0]
    outname = outpath + '/' + inname + '_paratope_segment.csv'
    segment_descdf.to_csv(outname, index=False)
    return outname

def plot_len_gaps():
    '''
    Plots gap properties
    '''
    infile = 'abdb_outfiles/respairs_gaps_cdr_numgaps.csv'
    df = pd.read_csv(infile)
    chains = ['H', 'L']
    f,axes = plt.subplots(2,1)
    for i,chain in enumerate(chains):
        ax  = axes[i]
        dfc = df[df.chain == chain]
        #sns.violinplot(x='cdr', y = 'plen', data=dfc, linewidth=0.0, ax=ax)
        sns.boxplot(x='cdr', y = 'plen', data=dfc, linewidth=0.2, ax=ax, fliersize=0)
        ax.set(xlabel = 'CDR %s'% chain)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    f.set_size_inches(7,7)
    plt.tight_layout()
    plotname = infile.split('/')[-1].split('.')[0] + '_plen'
    save_pdf(plotname) #plots paratopes length distribution
    plt.close()
    f, axes = plt.subplots(2,1)
    for i,chain in enumerate(chains):
        ax  = axes[i]
        dfc = df[df.chain == chain]
        #sns.violinplot(x='ngaps', y = 'plen', data=dfc, ax=ax, hue='cdr', linewidth=0)
        sns.boxplot(x='ngaps', y = 'plen', data=dfc, ax=ax, hue='cdr', linewidth=0.2, fliersize=0)
        ax.set(xlabel = '# of gaps %s'% chain)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    f.set_size_inches(7,7)
    plt.tight_layout()
    plotname = infile.split('/')[-1].split('.')[0] + '_plen_ngaps'
    save_pdf(plotname)

def plot_plen_sgaps():
    '''
    Plots length versus size of gaps
    '''
    infile = 'abdb_outfiles/respairs_gaps_cdr_numgaps.csv'
    df = pd.read_csv(infile)
    print(df.head(40))
    chains = ['H', 'L']
    f,axes = plt.subplots(2,1)
    for i,chain in enumerate(chains):
        ax  = axes[i]
        dfc = df[df.chain == chain]
        dfc= dfc[dfc.sgap>0]
        sns.boxplot(x='sgap', y='plen', data=dfc,ax=ax, hue='cdr', palette='pastel', linewidth=0.1, fliersize=0)
        #ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
        #ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        #ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
        ax.set(xlabel = 'sgap %s'% chain)
    f.set_size_inches(8,8)
    plt.tight_layout()
    plotname = infile.split('/')[-1].split('.')[0] + '_plen_sgap'
    save_pdf(plotname)
    plt.close()
    ## just CDR3
    dfc3 = df[df.cdr==3]
    dfc3 = dfc3[df.sgap > 0]
    sns.countplot(x='sgap', data=dfc3,linewidth=0.1)
    plotname = infile.split('/')[-1].split('.')[0] + '_sgap'
    save_pdf(plotname)

def plot_patterns():
    '''
    Investigates patterns emebeded in the paratopes (kmer thing)
    '''
    infile = 'abdb_outfiles/respairs_gaps_cdr_numgaps2.csv'
    df = pd.read_csv(infile)
    print(df.head())
    porder = df.shiftset.value_counts().index
    sns.set(style='dark')
    g = sns.FacetGrid(df, col='chain',row='cdr', margin_titles=True)
    g.map(sns.countplot, 'shiftset', order=df.shiftset.value_counts().iloc[:10].index, palette='pastel', linewidth=0)
    loc, labels= plt.xticks()
    g.set_xticklabels(labels=labels, rotation=90)
    plt.tight_layout()
    plotname = infile.split('/')[-1].split('.')[0] + '_shiftset'
    save_pdf(plotname)
    plt.close()
    g = sns.FacetGrid(df, col='chain',row='cdr', margin_titles=True)
    g.map(sns.countplot, 'gapset', order=df.gapset.value_counts().iloc[:10].index, palette='pastel', linewidth=0)
    loc, labels= plt.xticks()
    g.set_xticklabels(labels=labels, rotation=90)
    plt.tight_layout()
    plotname = infile.split('/')[-1].split('.')[0] + '_gapset'
    save_pdf(plotname)


def plot_patterns_segment():
    '''
    Investigates patterns emebeded in the paratopes (kmer thing)
    '''
    infile = 'abdb_outfiles/respairs_paratope_segment.csv'
    df = pd.read_csv(infile)
    print(df.head(40))
    print(df.shape)
    n = 25
    frac = df.gapset.value_counts()[:n].sum()/df.gapset.value_counts().sum()
    print(frac)
    print(len(df.gapset.unique()))
    topn  = df.gapset.value_counts().to_frame().iloc[:n]
    topn.columns = ['']
    #sns.heatmap(df.gapset.value_counts().to_frame().iloc[:n], yticklabels=1, cmap='Greys',
    sns.heatmap(topn, yticklabels=1, cmap='Greys',
    annot=True,fmt='d', linewidth=0)
    ax = plt.gca()
    ax.set_ylabel('gap patterns')
    plt.tight_layout()
    #view_temp()
    plotname = 'gap_patterns_top%s' % n
    save_pdf(plotname)
    plt.close()
    gapset_segment = pd.crosstab(df.gapset,df.segment)
    plt.figure(figsize=(15,55))
    sns.heatmap(gapset_segment, cmap='Greys', yticklabels=True, annot=True, fmt = 'd')
    ax = plt.gca()
    ax.set_ylabel('gap patterns')
    plt.tight_layout()
    plotname = 'gap_patterns_world'
    save_pdf(plotname)
    #return gapset_segment

def plot_shift_counts():
    '''
    Plots the shift counts
    '''
    infile = 'abdb_outfiles/respairs_gaps_cdr.csv'
    df = pd.read_csv(infile)
    print(df[df.pdbid == '4RGO_1'].iloc[40:60])
    sys.exit()
    order = df.abshift.value_counts().index[1:7]
    sns.set(style='dark')
    g = sns.FacetGrid(data=df, row='cdr', col='abchain', col_order=['H','L'],margin_titles=True)
    g.map(sns.countplot, 'abshift',palette='pastel',linewidth=0, order=order)
    plotname = infile.split('/')[-1].split('.')[0] + '_abshift'
    plt.tight_layout()
    save_pdf(plotname)


def get_pdb_infos():
    '''
    Grabs infos from the pdbfiles
    '''
    file_paths = get_single_antigens()
    labels = ['pdbid', 'lmol','lspec','hmol','hspec','agmol','agspec']
    contents  = []
    for fp in file_paths:
        pdbid = fp.split('/')[-1].split('.')[0]
        content = open(fp).read()
        parts = content.split('SEQRES')
        infos = [pdbid] + [item.split(':')[-1] for item in parts[0].splitlines()[-6:]]
        contents.append(infos)
    df = pd.DataFrame(contents, columns=labels)
    df.to_csv('abdb_outfiles/pdb_infos.csv', index=False)
    sys.exit()

def plot_species():
    '''
    Plots ab and ag species of the samples
    '''
    infile = 'abdb_outfiles/pdb_infos.csv'
    df = pd.read_csv(infile)
    vars = ['hspec', 'agspec']
    for var in vars:
        order = df[var].value_counts().index[:5]
        sns.countplot(df[var], order=order)
        loc, labels= plt.xticks()
        ax = plt.gca()
        ax.set_xticklabels(labels=labels, rotation=90)
        plotname = infile.split('/')[-1].split('.')[0] + '_' + var
        plt.tight_layout()
        save_pdf(plotname)
        plt.close()


def plot_plenset():
    '''
    Reformats respairs_gaps_cdr_numgaps.csv.
    '''
    infile = 'abdb_outfiles/respairs_gaps_cdr_numgaps2.csv'
    df = pd.read_csv(infile)
    plenset = []
    chains = ['H', 'L']
    vars = ['plen', 'plenresnum', 'maxgap', 'ngaps', 'cdr', 'chain']
    vals = []
    for name in df.name.unique():
        tempdf = df[df.name == name]
        for chain in chains:
            tempvals = [name]
            tempdf2 = tempdf[tempdf.chain == chain]
            try:
                for var in vars:
                    tempvals.append('-'.join([str(item) for item in tempdf2[var].tolist()]))
                tempvals[-1]  = tempvals[-1][0] # chain ID
                vals.append(tempvals)
            except:
                #print('attempt failed for: %s' % name)
                pass
    columns = ['pdbid'] + vars
    newdf = pd.DataFrame(vals, columns=columns)
    print(newdf[newdf.chain =='H'].plen.value_counts().shape)
    #g = sns.FacetGrid(data=newdf, col='chain')
    #g.map(sns.heatmap, 'plen')
    #view_temp()

def plot_abresnumiset(infile):
    '''
    check positions
    '''
    #infile = 'abdb_outfiles/respairs_paratope_segment.csv'
    df = pd.read_csv(infile)
    abresnums = []
    for item in df.abresnumiset:
        parts = item.split('-')
        abresnums += parts
    abresnumsdf = pd.DataFrame(abresnums, columns=['position'])
    plt.figure(figsize=(20,5))
    n = 50
    order = natsort.natsorted(abresnumsdf.position.value_counts().iloc[:n].index.sort_values())
    sns.countplot(abresnumsdf.position, order =order)
    loc, labels= plt.xticks()
    ax = plt.gca()
    ax.set_xticklabels(labels=labels, rotation=90)
    plt.tight_layout()
    plotname = 'position_distribution_top%s' % n
    save_pdf(plotname)
    plt.close()
    abresnums_segment = []
    for i, row in df.iterrows():
        abresnums = row.abresnumiset.split('-')
        for abresnum in abresnums:
            abresnums_segment.append((row.segment ,abresnum))
    abresnumsdf = pd.DataFrame(abresnums_segment, columns = ['segment', 'position'])
    abresnums_segment = pd.crosstab(index = abresnumsdf.position, columns = abresnumsdf.segment)
    print(abresnums_segment.head(50))
    print(abresnums_segment.info())
    order = natsort.natsorted(abresnums_segment.index)
    print(order)
    #abresnums_segment = abresnums_segment[order]
    abresnums_segment = abresnums_segment.reindex(order)
    plt.figure(figsize=(10,18))
    sns.heatmap(abresnums_segment, cmap='Greys', yticklabels=1, annot=True, fmt='d')
    plt.tight_layout()
    save_pdf('frequency_position')

def plot_residue_position_cdr3(infile):
    '''
    checks some stuff
    '''
    #abshiftloop =  'abdb_outfiles/respairs_absort_abresnumi_segments_abshift_abshiftl.csv'
    #respairs_segment = 'abdb_outfiles/respairs_paratope_segment.csv'
    df = pd.read_csv(infile)
    pd.set_option('display.max_column', None)
    df.iloc[:10].to_latex('abdb_tables/respairs_paratope_segment_snippet.tex')
    df.iloc[:10].to_csv('abdb_tables/respairs_paratope_segment_snippet.csv', index=False)
    sys.exit()
    cdr3df = df[df.segment == 'CDR-H3']
    # cdr3df = df
    paratope_resnums = []
    for i, row in cdr3df.iterrows():
        if str(row.paratope) != str(np.nan):
            print(row.paratope)
            paratope_residues = list(row.paratope)
            resnums = row.abresnumiset.split('-')
            paratope_resnum = [(residue, resnum) for residue, resnum in zip(paratope_residues, resnums)]
            paratope_resnums += paratope_resnum
    print(df.head())

    prdf = pd.DataFrame(paratope_resnums, columns= ['residue', 'position'])
    prcrossdf = pd.crosstab(prdf.position, prdf.residue)
    sindex = natsort.natsorted(prcrossdf.index)
    prcrossdf = prcrossdf.reindex(sindex)
    plt.figure(figsize=(10,10))
    sns.heatmap(prcrossdf, annot=True, fmt='d', cmap='Greys')
    outname= 'residue_position_cdr3'
    view_temp()
    sys.exit()
    save_pdf(outname)

def get_abseq_full(infiles):	
    '''
    Gets anibody sequence along with residue positions
    input:
    ------
    infiles = path to input files
    '''
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
    df.to_csv('abdb_outfiles/abdb_segment.csv',index=False)


def prepdata_ab():
    '''
    prepares data tables for antibodies.
    '''
    single_antigens = get_single_antigens() #get pdb with single antigen
    absorted, agsorted = get_residue_pairs_ab2(single_antigens[:]) # sort by antibody
    abinsert = get_unique_abresnumi(absorted) #account for inserted residues
    absegment = add_segments(abinsert) #add segments based on Martin numbering
    abshift = add_abshift(absegment) # add shift
    abshiftloop = add_abshiftl(abshift) #add shifft loop wise
    gap_patterns = get_numgaps_segment(abshiftloop)	#get gap patterns data
    #get_abseq_full(single_antigens)


def add_agshift(infile, outpath):
    '''
    adds shift sets from the antigen side
    :return:
    '''
    pd.set_option('display.max_column',None)
    df = pd.read_csv(infile)
    pdbids = df.pdbid.unique()
    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    # ids, abchains, segs, epis, elens, shifts, gaps, agres = [], [], [], [], [], [], [], []
    outdata = []
    excepts = []
    for pdbid in pdbids:
        # try:
        pdbdf = df[df.pdbid == pdbid]
        segments = pdbdf.segment.unique()
        for segment in segments:
            segdf = pdbdf[pdbdf.segment == segment]
            segdf = segdf.drop_duplicates(subset='agresnumi')
            #disregard residue order
            # sorted_residues = sorted(np.array(segdf.agresnumi.unique(), dtype=int))
            sorted_residues = np.array(segdf.agresnum.unique(), dtype=int)
            print(sorted_residues)
            # if pdbid == '4YPG_1':
            #     print(segdf)
            # calculate absolute distance (neglect ag residue order)
            if len(sorted_residues) > 1:
                shiftset = '-'.join([str(abs(sorted_residues[i]-sorted_residues[i-1])) for i in range(1,len(sorted_residues))])
                gapset = '-'.join([str(abs((sorted_residues[i]-sorted_residues[i-1]))-1) for i in range(1,
                                    len(sorted_residues))])
                agresnumiset = '-'.join([str(abs(item)) for item in sorted_residues])
            else:
                shiftset = np.nan
                gapset = np.nan
                agresnumiset = '-'.join([str(abs(item)) for item in sorted_residues])
            # print(shiftset)
            # print(gapset)
            abchain = segdf.abchain.unique()[0]
            agchain = segdf.agchain.unique()[0]
            # epitope = ''.join(aadict[aa] for aa in segdf.agres.values)
            segdf = segdf.sort_values(by = 'agresnum')
            epitope = ''.join(aadict[aa] for aa in segdf.agres.values)
            elen = len(epitope)
            outdata.append([pdbid, abchain, segment, epitope, elen, shiftset, gapset, agresnumiset, agchain])
    #     except:
    #         excepts.append(pdbid)
    #         print(pdbid)
    # print(len(excepts), 'exceptedids')
    colnames = ['pdbid', 'abchain','segment','epitope','elen','eshiftset','egapset','agresnumiset', 'agchain']
    # outdata = [ids, abchains, segs, epis, elens, shifts, gaps, agres]
    outdf = pd.DataFrame(outdata, columns=colnames)
    print(outdf.head(30))
    print(len(outdf.pdbid.unique()))
    inname = infile.split('/')[-1].split('.')[0]
    outname = outpath + '/' + inname + '_epitope_segment.csv'
    # outfile = 'abdb_outfiles/respairs_epitope_segment.csv'
    outdf.to_csv(outname, index=False)

def prepdata_ag():
    '''
    prepares data tables for antibodies.
    '''
    add_agshift('abdb_outfiles/respairs_absort_abresnumi_segments_abshift_abshiftl.csv')



def plotdata():
    '''
    Plots the data
    '''
    #plot_abresnumiset('abdb_outfiles/respairs_paratope_segment.csv')
    #plot_residue_position_cdr3('abdb_outfiles/respairs_paratope_segment.csv')
    plot_patterns_segment()

def check_dewitt():
    '''
    checks dewitt data set for multiple Cs
    '''
    infile = '../datasets/dewitt_naive_cdr3_aa_15.csv'
    df = pd.read_csv(infile).loc[:1000,:]
    aais = []
    seqs = []
    for i, row in df.iterrows():
        aas = list(row.cdr3_aa_15)
        c_count = 0
        for pair in enumerate(aas):
            aais.append(pair)
            if pair[1] == 'C':
                c_count += 1
        if c_count > 1:
            seqs.append(row.cdr3_aa_15)
    seqsdf = pd.DataFrame(seqs, columns=['cdr3_aa_15'])
    seqsdf.to_csv('../datasets/dewitt_naive_cdr3_aa_15_multic.csv',index=False)
    df2 = pd.DataFrame(aais, columns=['position','residue'])
    position_residue = pd.crosstab(df2.position, df2.residue)
    plt.figure(figsize=(10,10))
    sns.heatmap(position_residue, annot=True, fmt='d', cmap='Greys')
    view_temp()

def check_abdb():
    '''
    looks into multiple Cs on abdb
    '''
    infile = 'abdb_outfiles/respairs_paratope_segment.csv'
    df = pd.read_csv(infile).dropna(subset=['paratope'])
    #df = df[df.paratope != np.nan]
    df = df[df.segment == 'CDR-H3']
    data = []
    seq_resnums = []
    for i, row in df.iterrows():
        seq = row.paratope
        resnums = row.abresnumiset
        #print(type(seq))
        seqs = list(seq)
        c_count = 0
        for aa in seqs:
            if aa == 'C':
                c_count += 1
        if c_count > 1:
            seq_resnums.append((row.pdbid,seq,resnums))
    seq_resnumsdf = pd.DataFrame(seq_resnums, columns=['pdbid','interacting residue', 'position'])
    #seq_resnumsdf.to_csv('abdb_outfiles/abbd_CDR-H3_multic.csv', index=False)
    abdb_segment_file = 'abdb_outfiles/abdb_segment.csv'
    abdb_segment = pd.read_csv(abdb_segment_file)
    print(abdb_segment.head())
    multic_pdbid = seq_resnumsdf.pdbid.tolist()
    print(multic_pdbid)
    #multic_fullseq = abdb_segment[abdb_segment['pdbid'].isin(multic_pdbid)]
    multic_fullseq = abdb_segment
    multic_frcdr3 = multic_fullseq[multic_fullseq.segment.isin(['CDR-H3','HFR3', 'HFR4'])]
    #multic_cdr3 = multic_fullseq[multic_fullseq.segment == 'HFR3']
    multic_frcdr3 = multic_frcdr3[multic_frcdr3['pdbid'].isin(multic_pdbid)]
    for pdbid in multic_frcdr3.pdbid.unique():
        pdbdf = multic_frcdr3[multic_frcdr3.pdbid == pdbid]
        for segment in pdbdf.segment.unique():
            segmentdf = pdbdf[pdbdf.segment==segment]
            seq = ''.join(segmentdf.aa_single)
            print(seq)
            resnum = '-'.join(segmentdf.resnumi)
            datum = (pdbid, segment, seq,resnum)
            data.append(datum)
    columns = ['pdbid', 'segment', 'sequence', 'resnum']
    df = pd.DataFrame(data, columns=columns)
    print(df.head(40))
    df.to_csv('abdb_outfiles/abdb_multic_frcdr3.csv')


def translate_pattern_cdr3():
    '''
    Translates gap patterns to language for CDR-H3
    '''
    paratope_segment = 'abdb_outfiles/respairs_paratope_segment.csv'
    paratope_segment = pd.read_csv(paratope_segment)
    shiftl = 'abdb_outfiles/respairs_absort_abresnumi_segments_abshift_abshiftl.csv'
    shiftl = pd.read_csv(shiftl)
    paratope_segment = paratope_segment[paratope_segment.segment == 'CDR-H3']
    res_num = []
    for i,row in paratope_segment.iterrows():
            seqs = list(row.paratope)
            resnums = row.abresnumiset.split('-')
            for pair in zip(seqs,resnums):
                res_num.append(pair)
    res_num = pd.DataFrame(res_num, columns=['residue', 'position'])
    order = natsort.natsorted(res_num.position.unique())
    res_num_cross = pd.crosstab(res_num.position, res_num.residue)
    res_num_cross = res_num_cross.reindex(order)
    #position_resfreq = {}
    counter = 0
    position_resfreq = []
    for position, row in res_num_cross.iterrows():
        resfreq = ','.join([str(freq)+res for freq,res in zip(row, row.index)])
        resfreq = '{%s}'%resfreq
        position_resfreq.append((counter, position, resfreq))
        counter += 1
    position_resfreq = pd.DataFrame(position_resfreq, columns = ['idx', 'position','resfreq'])
    position_resfreq.to_csv('abdb_outfiles/position_resfreq.csv', index=False)
    print(position_resfreq.head(30))
    print(paratope_segment.head())

def check_data():
    '''
    Perfroms some checks on datasets
    '''
    #check_dewitt()
    #check_abdb()
    #translate_pattern_cdr3()


def plot_residue_number_distribution(infile, outfile):
    '''
    plots the residue number distribution.
    '''
    df = pd.read_csv(infile)
    print(df.info())
    sns.set(palette='pastel', style='darkgrid')
    sns.violinplot(df.abresnum,df.abchain, linewidth=0.0)
    ax = plt.gca()
    ax.set_ylabel('Chain')
    ax.set_xlabel('Residue number')
    #plotname = infile.split('/')[-1].split('.')[0] + '_abresnum_distribution'
    save_pdf(outfile)


def plot_interacting_residue_distribution(infile, outfile):
    '''
    plots residue distributions and gaps
    infile: respair_gaps.csv
    '''
    df = pd.read_csv(infile)
    dfl = df[df.abchain == 'L'].abres.value_counts()
    dfh = df[df.abchain == 'H'].abres.value_counts()
    aagroups = [['ASP', 'GLU', 'LYS', 'ARG','SER', 'THR', 'GLN', 'ASN', 'HIS', 'CYS', 'ALA', 'VAL','LEU', 'ILE', 'PHE', 'TRP', 'MET','PRO','GLY','TYR'],
                ['c','c', 'c', 'c','p', 'p','p','p','p','p','h','h','h','h','h','h','h','h', 'h', 'p']]
    aagroups = pd.DataFrame(np.array(aagroups).T, columns=['aaname','group'])
    aagroups = aagroups.set_index('aaname')
    dfl = pd.concat([dfl.to_frame(),aagroups],axis=1)
    dfh = pd.concat([dfh.to_frame(),aagroups],axis=1)
    dfs = [dfh,dfl]
    f, axes = plt.subplots(1,len(dfs))
    chains = ['H','L']
    for i, chain in enumerate(chains):
        ax = axes[i]
        df = dfs[i]
        sns.barplot(y=df.index, x=df.abres,hue=df.group,palette='Greys',dodge=False, ax=ax)
        ax.set_xlabel('Counts in %s' % chain)
        ax.set(xlim=[0,2000])
    plt.tight_layout()
    save_pdf(outfile)


def plot_interacting_residue_distribution_surface(infile, outfile):
    '''
    plots interacting residues normalized by surface residues
    '''
    df = pd.read_csv(infile)
    dfl = df[df.abchain == 'L'].abres.value_counts()
    dfh = df[df.abchain == 'H'].abres.value_counts()
    surfacel = pd.read_csv('../datasets/surface_residues/surface_residue_count_L.csv')
    surfaceh = pd.read_csv('../datasets/surface_residues/surface_residue_count_H.csv')
    surfaces = [surfaceh, surfacel]
    aagroups = [['ASP', 'GLU', 'LYS', 'ARG','SER', 'THR', 'GLN', 'ASN', 'HIS', 'CYS', 'ALA', 'VAL','LEU', 'ILE', 'PHE', 'TRP', 'MET','PRO','GLY','TYR'],
                ['c','c', 'c', 'c','p', 'p','p','p','p','p','h','h','h','h','h','h','h','h', 'h', 'p']]
    aagroups = pd.DataFrame(np.array(aagroups).T, columns=['aaname','group'])
    aagroups = aagroups.set_index('aaname')
    dfl = pd.concat([dfl.to_frame(),aagroups],axis=1)
    dfh = pd.concat([dfh.to_frame(),aagroups],axis=1)
    dfs = [dfh,dfl]
    f, axes = plt.subplots(1,len(dfs))
    chains = ['H','L']
    for i, chain in enumerate(chains):
        ax = axes[i]
        df = dfs[i]
        df['normalized_count'] = (df['abres']-df.abres.min())/(df.abres.max()-df.abres.min())
        sdf = surfaces[i].set_index(surfaces[i].residue)
        print(sdf.head())
        sdf['normalized_count'] = (sdf.rescount-sdf.rescount.min())/(sdf.rescount.max()-sdf.rescount.min())
        df['normalized_surface_ratio'] = df.normalized_count/sdf.normalized_count
        df['normalized_surface'] = sdf.normalized_count
        #pd.set_option('display.max_column',None)
        #print(df)
        sns.barplot(y=df.index, x=df.normalized_surface_ratio,hue=df.group,palette='Greys',dodge=False, ax=ax)
        ax.set_xlabel('Interacting/surface %s' % chain)
        ax.set(xlim=[0,5])
    plt.tight_layout()
    save_pdf(outfile)

def plot_paratope_length(infile,outfile):
    '''
    plots the length of paratopes for each segment
    :return:
    '''
    sns.set(style='dark')
    df = pd.read_csv(infile)
    pd.set_option('display.max_column', None)
    print(df.head())
    plt.figure(figsize=(12,8))
    medians = df.groupby(['segment'])['plen'].agg(np.median)
    median_labels = [str(median) for median in medians]
    # ax = sns.boxplot(df.segment,df.plen, linewidth=0.1,fliersize=0, color='black',order=medians.index)
    ax = sns.boxplot(df.segment,df.plen, fliersize=0, color='lightgrey',order=medians.index)
    # ax = sns.boxplot(df.segment,df.plen, linewidth=0.1,fliersize=0, palette='Greys',order=medians.index)
    ax.set_ylabel('Paratope size')
    ax.set_xlabel('Segment')
    ax.set_ylim(0,12)
    pos = range(len(medians))
    pos_adj = [0,12]
    for tick, label in zip(pos, ax.get_xticklabels()):
        ax.text(pos[tick]+0.2, medians[tick]+0.1, median_labels[tick],
                horizontalalignment='center', size='x-large', color='grey', weight='bold')
        # if tick not in pos_adj:
        # 	ax.text(pos[tick], medians[tick], median_labels[tick],
        #             horizontalalignment='center', size='x-small', color='white', weight='bold')
        # else:
        # 	ax.text(pos[tick], medians[tick]-0.5, median_labels[tick],
        # 			horizontalalignment='center', size='x-small', color='white', weight='bold')

    save_pdf(outfile)


def plot_paratope_gapsize(infile,outfile):
    '''
    Plots number of gaps in paratopes, segment wise
    :param infile:
    :param outfile:
    :return:
    '''
    pd.set_option('display.max_column',None)
    df = pd.read_csv(infile)
    segment_gaps = []
    for i,row in df.iterrows():
        segment = row.segment
        try: # handles the nan in the dataset
            gaps = row.gapset.split('-')
            segment_gap = [(segment,int(gap)) for gap in gaps]
            for sg in segment_gap:
                segment_gaps.append(sg)
        except AttributeError:
            print('found a nan, skipping...')
    sgdf = pd.DataFrame(segment_gaps, columns=['Segment','Gap size'])
    plt.figure(figsize=(12,8))
    sns.set(style='dark', palette='Greys')
    order = sorted(sgdf.Segment.unique())
    segment_max = []
    for ord in order:
        orddf = sgdf[sgdf.Segment == ord]
        maxlen = orddf.max()['Gap size']
        segment_max.append((ord,maxlen))
    segment_max = pd.DataFrame(segment_max, columns=['Segment', 'Maximum gap'])
    segment_max.to_latex('../tex/tables/max_gap.tex', index=False)
    sys.exit()
    sns.stripplot(x='Segment', y='Gap size',data=sgdf, color='black', order=order,alpha=0.2, linewidth=0, size=13)
    save_pdf(outfile)

def plot_gap_patterns_gapset(infile,outfile):
    '''
    Plots gap patterns found on gapset
    :param infile:
    :param outfile:
    :return:
    '''
    pd.set_option('display.max_column', None)
    df = pd.read_csv(infile)
    crosstabdf = pd.crosstab(index=df.gapset, columns=df.segment)
    crosstabdf_norm  = (crosstabdf - crosstabdf.min())/(crosstabdf.max()-crosstabdf.min())
    crosstabdf2 = crosstabdf_norm.loc[(crosstabdf_norm > 0.2).any(1)]
    linear_patterns = [item for item in crosstabdf2.index if '0-0' in item]
    fraction_linear_patterns = len(linear_patterns)/len(crosstabdf2.index)
    print(fraction_linear_patterns)
    sys.exit()
    #plt.figure(figsize=(15,70))
    plt.figure(figsize=(15,15))
    sns.heatmap(crosstabdf2, cmap='Greys', yticklabels=True,annot=True)
    ax = plt.gca()
    ax.set_ylabel('Gap patterns')
    ax.set_xlabel('Segment')
    plt.tight_layout()
    save_pdf(outfile)

def add_gap_patterns_position_residue(infile,outfile):
    '''

    :param infile: a csv file
    :param outfile: a csv file
    :return:
    '''
    pd.set_option('display.max_columns', None)
    df = pd.read_csv(infile).dropna()
    gap_position_set = []
    gap_residue_position_set = []
    for i,row in df.iterrows():
        gaps = row.gapset.split('-')
        positions = row.abresnumiset.split('-')
        residues = list(row.paratope)
        gap_position = ''
        gap_residue_position = ''
        for i2,gap in enumerate(gaps):
            gap_position += positions[i2]
            gap_position += '[%s]' % gap
            gap_residue_position += residues[i2]
            gap_residue_position += positions[i2]
            gap_residue_position += '[%s]' % gap
        gap_position += positions[i2+1]
        gap_residue_position += residues[i2+1]
        gap_residue_position += positions[i2+1]
        gap_position_set.append(gap_position)
        gap_residue_position_set.append(gap_residue_position)
    df['gap_position_set'] = gap_position_set
    df['gap_residue_position_set'] = gap_residue_position_set
    df.to_csv(outfile, index=False)


def plot_gap_patterns_position(infile,outfile):
    '''
    Plots gap patterns found on gapset
    :param infile:
    :param outfile:
    :return:
    '''
    pd.set_option('display.max_column', None)
    df = pd.read_csv(infile)
    # df = df[df.abchain=='H']
    crosstabdf = pd.crosstab(index=df.gapset, columns=df.segment)
    crosstabdf_norm  = (crosstabdf - crosstabdf.min())/(crosstabdf.max()-crosstabdf.min())
    crosstabdf2 = crosstabdf_norm.loc[(crosstabdf_norm > 0.2).any(1)]
    gap_patterns = crosstabdf2.index
    all_residues = []
    all_positions = []
    for gap_pattern in gap_patterns[:]:
        residues = []
        positions = []
        gapdf = df[df.gapset==gap_pattern]
        for i,row in gapdf.iterrows():
            paratope_residues = list(row.paratope)
            for residue in paratope_residues:
                residues.append(residue)
                positions.append(row.abresnumiset)
        all_positions += positions
        all_residues += residues
    resposdf = pd.DataFrame({'Residue':all_residues,'Position':all_positions})
    ordered_index = natsort.natsorted(resposdf.Position.unique())
    resposct = pd.crosstab(resposdf.Position, resposdf.Residue)
    resposct = resposct.reindex(ordered_index)
    resposct = (resposct-resposct.min())/(resposct.max()-resposct.min())
    resposct = resposct.round(1)
    resposct = resposct.loc[(resposct > 0.4).any(1)]
    # resposct = resposct.loc[(resposct > 30).any(1)]
    #plt.figure(figsize=(15,70))
    print(len(gap_patterns))
    print(len(resposct.index))
    plt.figure(figsize=(14,10))
    sns.heatmap(resposct, cmap='Greys', yticklabels=True,annot=True, square=False)
    ax = plt.gca()
    ax.set_ylabel('Position')
    ax.set_xlabel('Residue')
    plt.tight_layout()
    # view_temp()
    # sys.exit()
    # outname = outfile + '_' + gap_pattern
    save_pdf(outfile)


def plot_gap_patterns_position_chain(infile,outfile):
    '''
    Plots gap patterns found on gapset
    :param infile:
    :param outfile:
    :return:
    '''
    pd.set_option('display.max_column', None)
    sdf = pd.read_csv(infile)
    # df = df[df.abchain=='H']
    chains = ['L','H']
    for chain in chains:
        df = sdf[sdf.abchain==chain]
        crosstabdf = pd.crosstab(index=df.gapset, columns=df.segment)
        crosstabdf_norm  = (crosstabdf - crosstabdf.min())/(crosstabdf.max()-crosstabdf.min())
        crosstabdf2 = crosstabdf_norm.loc[(crosstabdf_norm > 0.2).any(1)]
        gap_patterns = crosstabdf2.index
        print(gap_patterns)
        all_residues = []
        all_positions = []
        for gap_pattern in gap_patterns[:]:
            residues = []
            positions = []
            gapdf = df[df.gapset==gap_pattern]
            for i,row in gapdf.iterrows():
                paratope_residues = list(row.paratope)
                for residue in paratope_residues:
                    residues.append(residue)
                    positions.append(row.abresnumiset)
            all_positions += positions
            all_residues += residues
        resposdf = pd.DataFrame({'Residue':all_residues,'Position':all_positions})
        ordered_index = natsort.natsorted(resposdf.Position.unique())
        resposct = pd.crosstab(resposdf.Position, resposdf.Residue)
        resposct = resposct.reindex(ordered_index)
        resposct = (resposct-resposct.min())/(resposct.max()-resposct.min())
        resposct = resposct.loc[(resposct > 0.3).any(1)]
        #plt.figure(figsize=(15,70))
        plt.figure(figsize=(15,15))
        sns.heatmap(resposct, cmap='Greys', yticklabels=True,annot=True)
        ax = plt.gca()
        ax.set_ylabel('Position')
        ax.set_xlabel('Residue')
        plt.tight_layout()
        outname = outfile + '_' + chain
        save_pdf(outname)
        plt.close()

def plot_gap_patterns_distribution(infile):
    '''
    Plots gap patterns distribution
    :return:
    '''
    df = pd.read_csv(infile)
    df = df[df.segment == 'CDR-H3']
    normalized_counts = df.gapset.value_counts()/sum(df.gapset.value_counts())*100
    n = 12
    top_n = normalized_counts.sort_values().iloc[-n:]
    top_6  = sum(top_n[-6:])
    print('top 6: %s' % top_6 )
    sns.barplot(top_n.index, top_n, color='lightgrey')
    plt.xticks(rotation=90)
    plt.tight_layout()
    outname1 = 'gap_pattern_percent_cdr3'
    plt.ylabel('percent (%)')
    save_pdf(outname1)
    plt.close()
    sys.exit()
    print(top_n.index)
    f, axes = plt.subplots(figsize=(25,20), nrows=4, ncols=3)
    axes = axes.flatten()
    for axi, gap_pattern in enumerate(top_n.index):
        gapdf = df[df.gapset == gap_pattern]
        paratope_resnums = []
        for i, row in gapdf.iterrows():
            if str(row.paratope) != str(np.nan):
                print(row.paratope)
                paratope_residues = list(row.paratope)
                resnums = row.abresnumiset.split('-')
                paratope_resnum = [(residue, resnum) for residue, resnum in zip(paratope_residues, resnums)]
                paratope_resnums += paratope_resnum
        prdf = pd.DataFrame(paratope_resnums, columns= ['residue', 'position'])
        prcrossdf = pd.crosstab(prdf.position, prdf.residue)
        sindex = natsort.natsorted(prcrossdf.index)
        prcrossdf = prcrossdf.reindex(sindex)
        ax = axes[axi]
        sns.heatmap(prcrossdf, annot=True, fmt='d', cmap='Greys', ax=ax)
        ax.set_xlabel(gap_pattern)
    outname2= 'gap_pattern_decomposition_cdr3'
    save_pdf(outname2)


def plot_gap_pos_patterns_distribution(infile):
    '''
    Plots gap patterns distribution
    :return:
    '''
    df = pd.read_csv(infile)
    df = df[df.segment == 'CDR-H3']
    print(df.head())
    normalized_counts = df.abresnumiset.value_counts()/sum(df.abresnumiset.value_counts())*100
    n = 12
    top_n = normalized_counts.sort_values().iloc[-n:]
    sns.barplot(top_n.index, top_n, color='lightgrey')
    plt.xticks(rotation=90)
    plt.tight_layout()
    outname1 = 'gap_pospattern_percent_cdr3'
    plt.ylabel('percent (%)')
    save_pdf(outname1)
    plt.close()
    print(top_n.index)
    f, axes = plt.subplots(figsize=(25,20), nrows=4, ncols=3)
    axes = axes.flatten()
    for axi, gap_pattern in enumerate(top_n.index):
        gapdf = df[df.abresnumiset == gap_pattern]
        paratope_resnums = []
        for i, row in gapdf.iterrows():
            if str(row.paratope) != str(np.nan):
                print(row.paratope)
                paratope_residues = list(row.paratope)
                resnums = row.abresnumiset.split('-')
                paratope_resnum = [(residue, resnum) for residue, resnum in zip(paratope_residues, resnums)]
                paratope_resnums += paratope_resnum
        prdf = pd.DataFrame(paratope_resnums, columns= ['residue', 'position'])
        prcrossdf = pd.crosstab(prdf.position, prdf.residue)
        sindex = natsort.natsorted(prcrossdf.index)
        prcrossdf = prcrossdf.reindex(sindex)
        ax = axes[axi]
        sns.heatmap(prcrossdf, annot=True, fmt='d', cmap='Greys', ax=ax)
        ax.set_xlabel(gap_pattern)
    outname2= 'gap_pospattern_decomposition_cdr3'
    save_pdf(outname2)


def pdb2seqtable():
    '''
    Transforms some pdbfiles to sequence and residue number table
    :return:
    '''
    file_paths = fifi.find_files('../datasets/NR_LH_Protein_Martin','.pdb')
    n = 5
    # fps = random.sample(file_paths,5)
    fps = file_paths
    data = {}
    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    print(aadict)
    for fp in fps:
        parser = PDB.PDBParser(PERMISSIVE=1)
        structure_name = fp.split('/')[-1].split('.')[0]
        structure = parser.get_structure(structure_name, fp)
        abchains = ['H', 'L']
        chains = PDB.Selection.unfold_entities(structure,'C')
        chains = [item for item in chains if item.id in abchains]
        if len(chains[0].get_list()) == len(chains[1].get_list()) and len(chains[0].get_list()) == 115:
            print(structure_name)
            for chain in chains:
                header = structure.id + '_' + chain.id
                contents = []
                if chain.id in abchains:
                    residues = chain.get_residues()
                    for residue in residues:
                        resname = residue.resname
                        resname_one = aadict[resname]
                        resseq = residue.full_id[-1][1]
                        contents.append('%s-%s-%s' % (resname_one,resname,resseq))
                data[header] = contents
    df = pd.DataFrame(data,)
    outfile = 'abdb_tables/martin_numbered_proteins.tex'
    df.to_latex(outfile, index=False)

def plot_gap_patterns_all_segments(infile):
    '''
    explores plotnine for future plottings
    :return:
    '''
    # pd.set_option('display.max_column', None)
    df = pd.read_csv(infile)
    n = -1
    topsdfs = []
    for segment in df.segment.unique():
        sdf = df[df.segment == segment]
        value_counts = sdf.gap_patterns.value_counts().iloc[:]
        top_patterns = value_counts.index
        total_pattern = sum(sdf.gap_patterns.value_counts())
        topsdf = sdf[sdf.gap_patterns.isin(top_patterns)]
        value_dict = {item:value for item,value in zip(value_counts.index, value_counts.values)}
        contributions = []
        for i,row in topsdf.iterrows():
            contribution = value_dict[row.gap_patterns]
            contributions.append(contribution)
            # print(contribution)
        topsdf['raw_count'] = contributions
        topsdf = topsdf.drop_duplicates(subset='gap_patterns')
        sorted_topsdf = topsdf.sort_values('raw_count', ascending=False)
        ranks = [i for i in range(1,topsdf.shape[0]+1)]
        sorted_topsdf['order'] = ranks
        sorted_topsdf['normalized_count'] = sorted_topsdf.raw_count/total_pattern
        gapstrstatus = []
        for status in sorted_topsdf.gapstatus:
            if status == 0:
                gapstrstatus.append('continuous')
            else:
                gapstrstatus.append('discontinuous')
        sorted_topsdf['gapstrstatus'] = gapstrstatus
        topsdfs.append(sorted_topsdf)
        # sizes = sorted_topsdf.raw_count.groupby(sorted_topsdf.gapstrstatus).sum()
        # plt.pie(sizes, autopct='%1.1f%%', labels=['continuous', 'discontinuous'])
        outpiename = infile.split('/')[1].split('.')[0] + '_top%s_gap_patterns_%s_pie.pdf'%(n,segment)
        # save_pdf(outpiename)
        # plt.close()
    topsdf_conc = pd.concat(topsdfs)
    outcsvname = 'abdb_outfiles/' + '_'.join(outpiename.split('.')[0].split('_')[:-2]) + '.csv'
    topsdf_conc.to_csv(outcsvname, index=False)
    print(outcsvname)
    sys.exit()
    gapstatus_cat = CategoricalDtype(categories=topsdf_conc['gapstrstatus'].unique().astype(str).tolist())
    topsdf_conc['gapstatus_category'] = topsdf_conc['gapstrstatus'].astype(str).astype(gapstatus_cat)
    gplot = (ggplot(topsdf_conc)
             # + aes('contribution', fill= 'gapstatus_ord')
             + aes('order', 'raw_count', fill='gapstatus_category')
             + geom_bar(stat='identity')
             # + geom_jitter(color='gapstatus_category')
             # + facet_grid(('absegment','abchain'), scales='free')
             + facet_wrap(('segment'), scales='free', ncol=3)
             + scale_fill_cmap_d('viridis')
             # + scale_x_discrete(breaks=range(1,800), labels = [str(i) for i in range(1,800)], limits=range(1,800))
             # + scale_x_discrete(breaks=range(1,500)
             # + scale_x_continuous(breaks='order', labels='gap_patterns',limits=(1,10))
             # + coord_flip()
             # + xlab('rank')
             + ylab('fraction in top %s'%n)
             + geom_text(aes(label='gap_patterns'),  # new
                         size=2, va='center', nudge_x=0, nudge_y=2 )
             + theme(panel_spacing = 0.5)
             + theme(axis_text_x=element_text(rotation=0), figure_size=(25,25))
             + theme(strip_background=element_rect(fill='#2266FF', size=1.4, alpha=.30),
                     panel_background=element_rect(fill='#2266FF', alpha=.05))
             # + scale_y_log10()
             )\
            + coord_flip()
    outname = 'abdb_figures/' + infile.split('/')[1].split('.')[0] + '_top%s_gap_patterns_all_segments.pdf' % n
    # gview_temp(gplot)
    # sys.exit()
    gplot.save(outname)
    os.system('open %s' % outname)



def batch_add_gap_count_data(infiles):
    '''
    add gap counts data
    :return:
    '''
    # pd.set_option('display.max_column', None)
    for infile in infiles:
        df = pd.read_csv(infile)
        n = -1
        topsdfs = []
        for segment in df.segment.unique():
            sdf = df[df.segment == segment]
            value_counts = sdf.gap_patterns.value_counts().iloc[:]
            top_patterns = value_counts.index
            total_pattern = sum(sdf.gap_patterns.value_counts())
            topsdf = sdf[sdf.gap_patterns.isin(top_patterns)]
            value_dict = {item:value for item,value in zip(value_counts.index, value_counts.values)}
            contributions = []
            for i,row in topsdf.iterrows():
                contribution = value_dict[row.gap_patterns]
                contributions.append(contribution)
                # print(contribution)
            topsdf['raw_count'] = contributions
            topsdf = topsdf.drop_duplicates(subset='gap_patterns')
            sorted_topsdf = topsdf.sort_values('raw_count', ascending=False)
            ranks = [i for i in range(1,topsdf.shape[0]+1)]
            sorted_topsdf['order'] = ranks
            sorted_topsdf['normalized_count'] = sorted_topsdf.raw_count/total_pattern
            topsdfs.append(sorted_topsdf)
        topsdf_conc = pd.concat(topsdfs)
        # outcsvname = 'abdb_outfiles/' + '_'.join(outpiename.split('.')[0].split('_')[:-2]) + '.csv'
        outcsvname = infile.split('.')[0] + '_count.csv'
        print(outcsvname)
        topsdf_conc.to_csv(outcsvname, index=False)



def plot_gap_pos_patterns_all_segments(infile):
    '''
    explores plotnine for future plottings
    :return:
    '''
    # pd.set_option('display.max_column', None)
    df = pd.read_csv(infile)
    pd.set_option('display.max_column', None)
    n = 25
    top_n = df.abresnumiset.value_counts().iloc[:n].index
    segments = df.segment.unique()
    cdf = pd.DataFrame()
    for segment in segments:
        tempdf = df[df.segment == segment]
        top3_patterns = tempdf.abresnumiset.value_counts().iloc[:3].index
        top3df = tempdf[tempdf.abresnumiset.isin(top3_patterns)]
        cdf = pd.concat([cdf,top3df])
    topdf = df[df.abresnumiset.isin(top_n)]
    gplot = (ggplot(topdf)
             + aes('abresnumiset', fill='abresnumiset')
             + geom_bar()
             + facet_grid('absegment~abchain', scales='free')
             + scale_fill_cmap_d('viridis')
             + theme(axis_text_x=element_text(rotation=90), figure_size=(10,20))
             + theme(strip_background=element_rect(fill='#2266FF', size=1.4, alpha=.30),
                     panel_background=element_rect(fill='#2266FF', alpha=.05))
             )
    outname = 'abdb_figures/top%s_gap_pos_patterns_all_segments.pdf' % n
    gplot.save(outname)
    os.system('open %s' % outname)

def add_notationx(infile, gapsettype,resnumitype):
    '''
    adds a new gap pattern notation. 0-0 > XXX, 2-2 > X2X2X
    also adds cdr columns without chain annotation CDR-L1 to CDR1
    also adds gapstatus: 0 no gap, 1 with gap
    :return:
    '''
    # infile = 'abdb_outfiles/respairs_paratope_segment.csv'
    df = pd.read_csv(infile)
    gap_patterns = []
    newsegments = []
    gapstats = []
    for i,row in df.iterrows():
        newsegment = re.sub('[LH-]','',row.segment)
        newsegments.append(newsegment)
        gapstatflag = 0 # a flag for gap status
        if len(row[resnumitype]) > 1 and str(row[gapsettype]) != 'nan':
            # print(row[gapsettype])
            gaps = row[gapsettype].split('-')
            gap_pattern = 'X'
            for gap in gaps:
                if gap == '0':
                    gap_pattern += 'X'
                else:
                    gap_pattern += '%sX' % gap
                    gapstatflag = 1
            gap_patterns.append(gap_pattern)
        else:
            gap_patterns.append('X')
        if gapstatflag:
            gapstats.append(1)
        else:
            gapstats.append(0)
    df['gap_patterns'] = gap_patterns
    df['absegment'] = newsegments
    df['gapstatus'] = gapstats
    outfile = infile.split('.')[0] + '_notationx.csv'
    df.to_csv(outfile, index=False)


def batch_add_notationx(infiles):
    '''
    batch add a new gap pattern notation. 0-0 > XXX, 2-2 > X2X2X
    also adds cdr columns without chain annotation CDR-L1 to CDR1
    also adds gapstatus: 0 no gap, 1 with gap
    :return:
    '''
    # infile = 'abdb_outfiles/respairs_paratope_segment.csv'
    for infile in infiles:
        gapsettype = 'gapset'
        resnumitype = 'abresnumiset'
        if 'epitope' in infile:
            gapsettype = 'egapset'
            resnumitype = 'agresnumiset'
        df = pd.read_csv(infile)
        gap_patterns = []
        newsegments = []
        gapstats = []
        gapstrstatus = []
        for i,row in df.iterrows():
            newsegment = re.sub('[LH-]','',row.segment)
            newsegments.append(newsegment)
            gapstatflag = 0 # a flag for gap status
            if len(row[resnumitype]) > 1 and str(row[gapsettype]) != 'nan':
                # print(row[gapsettype])
                gaps = row[gapsettype].split('-')
                gap_pattern = 'X'
                for gap in gaps:
                    if gap == '0':
                        gap_pattern += 'X'
                    else:
                        gap_pattern += '%sX' % gap
                        gapstatflag = 1
                gap_patterns.append(gap_pattern)
            else:
                gap_patterns.append('X')
            if gapstatflag:
                gapstats.append(1)
                gapstrstatus.append('discontinuous')
            else:
                gapstats.append(0)
                gapstrstatus.append('continuous')
        df['gap_patterns'] = gap_patterns
        df['absegment'] = newsegments
        df['gapstatus'] = gapstats
        df['gapstrstatus'] = gapstrstatus
        outfile = infile.split('.')[0] + '_notationx.csv'
        df.to_csv(outfile, index=False)

def get_pattern_residue_decomposition(infile, parepi, resnumiset, absegment):
    '''
    plots the residue-position decompostions of the gap patterns
    :return:
    '''
    df = pd.read_csv(infile)
    if absegment != 'all':
        df = df[df.segment == absegment]
        face_wrap_ncol = 3
    else:
        face_wrap_ncol = 5
    n = 10
    excepts = 0
    topn_gap_pattern = df.gap_patterns.value_counts().iloc[:n].index.tolist()
    newdf = {'residue':[], 'position':[],'gap_pattern':[], 'frequency':[]}
    for gap_pattern in topn_gap_pattern:
        gpdf = df[df['gap_patterns'] == gap_pattern]
        residues = []
        positions = []
        for i, row in gpdf.iterrows():
            try:
                parepis = list(row[parepi])
                residues += parepis
                pos  = row[resnumiset].split('-')
                positions += pos
                if not len(parepis) == len(pos): print(row)
            except:
                excepts += 1
        ctdf = pd.DataFrame()
        ctdf['residue'] = residues
        ctdf['position'] = positions
        crossed = pd.crosstab(ctdf.position, ctdf.residue)
        aacids = crossed.columns
        nrows = crossed.shape[0]
        for aacid in aacids:
            frequency = crossed[aacid]
            newdf['frequency'] += frequency.values.tolist()
            newdf['position'] += frequency.index.values.tolist()
            newdf['gap_pattern'] += [gap_pattern]*nrows
            newdf['residue'] += [aacid]*nrows
    print('failed parepis: %s' % excepts)
    tempdf = pd.DataFrame(newdf)
    outdf = infile.split('.')[0] + '_respos_decomposition_%s.csv' % absegment
    tempdf.to_csv(outdf, index=False)
    # topn_pos = tempdf.position.value_counts().iloc[:n].index.tolist()
    # topn_pos = tempdf.position.value_counts().iloc[:n].index.tolist()
    # print(topn_pos, 'heyy im here')
    # topn_pos = natsort.natsorted(topn_pos)
    # topn_pos_cat = CategoricalDtype(categories=topn_pos, ordered=True)
    # tempdf = tempdf[tempdf.position.isin(topn_pos)]
    # tempdf['position_ord'] = tempdf['position'].astype(str).astype(topn_pos_cat)
    tempdf['log_frequency'] = np.log(tempdf.frequency)
    #drop zeros for visibility?
    tempdf = tempdf.sort_values(by='frequency').iloc[-100:]
    if 'epitope' in infile:
        tempdf['position'] = tempdf['position'].astype(int)
        g = ( ggplot(tempdf)
              + aes('residue','position', fill='frequency')
              + geom_tile(aes(width=0.99, height=0.99))
              + facet_wrap('gap_pattern', ncol=face_wrap_ncol, scales='free')
              + coord_equal()
              + geom_text(aes(label='frequency'), color='white')
              # + theme(strip_background=element_rect(fill='#2266FF', size=1.4, alpha=.30),
              #         panel_background=element_rect(fill='#440154FF', alpha=1), line=element_line(color='#440154FF'))
                      # panel_background=element_rect(fill='#400b0b', alpha=1), line=element_line(color='#400b0b'))
              + theme(figure_size=(20,25), panel_spacing=1)
              # + ylab('log position')
              # + scale_y_log10()
            )
    elif 'paratope' in infile:
        g = ( ggplot(tempdf)
              + aes('residue','position', fill='frequency')
              + geom_tile(aes(width=0.99, height=0.99))
              + facet_wrap('gap_pattern', ncol=face_wrap_ncol, scales='free')
              + coord_equal()
              + geom_text(aes(label='frequency'), color='white')
              # + theme(strip_background=element_rect(fill='#2266FF', size=1.4, alpha=.30),
              #         panel_background=element_rect(fill='#440154FF', alpha=1), line=element_line(color='#440154FF'))
              # panel_background=element_rect(fill='#400b0b', alpha=1), line=element_line(color='#400b0b'))
              + theme(figure_size=(20,25), panel_spacing=1)
              # + scale_y_log10()
              )
    # outfig = 'abdb_figures/gap_pattern_residue_position_decomposition.pdf'
    outfig = 'abdb_figures/' + infile.split('/')[-1].split('.')[0] + '_respos_decomposition_%s.pdf' % absegment
    g.save(outfig)
    os.system('open %s' % outfig)

def plot_epitope_length(infile):
    '''
    plots the distribution of epitope lengths
    :param infile:
    :return:
    '''
    df = pd.read_csv(infile)
    print(df.head())
    g = ( ggplot(df)
          + aes('absegment', 'elen')
          # + geom_boxplot()
          + geom_violin()
          + facet_wrap('abchain')
          + theme(figure_size=(7,6))
    )
    outfig = infile.split('/')[-1].split('.')[0] + '_elen_boxplot'
    # gview_temp(g)
    # sys.exit()
    gsave_pdf(g, outfig)



def plot_epitope_gapsize(infile):
    '''
    plots the distribution of epitope lengths
    :param infile:
    :return:
    '''
    df = pd.read_csv(infile)
    print(df.head())
    tempcols = []
    for i,row in df.iterrows():
        try:
            gapset = row.egapset.split('-')
            for gap in gapset:
                tempcols.append((row.abchain, row.absegment, int(gap)))
        except:
            print(' %s found, skipping this one' % str(row.egapset))
    tempdf = pd.DataFrame(tempcols, columns=['abchain', 'absegment', 'gapsize'])
    g = ( ggplot(tempdf)
          + aes('absegment', 'gapsize')
          + geom_violin()
          + facet_wrap('abchain')
          + theme(figure_size=(7,6))
          + scale_y_log10()
          + ylab('log gapsize')
          )
    outfig = infile.split('/')[-1].split('.')[0] + '_gapsize_boxplot'
    gsave_pdf(g, outfig)

def plot_pattern_residue_decomposition_cdrs_ab():
    '''
    plots respos decomposition for all cdr segments in ab
    :return:
    '''
    segments = ['CDR-H3', 'CDR-H2', 'CDR-H1', 'CDR-L3', 'CDR-L2', 'CDR-L1']
    for segment in segments:
        get_pattern_residue_decomposition('abdb_outfiles/respairs_paratope_segment_notationx.csv', 'paratope',
                                          'abresnumiset', segment)


def plot_pattern_residue_decomposition_cdrs_ag():
    '''
    plots respos decomposition for all cdr segments in ab
    :return:
    '''
    segments = ['CDR-H3', 'CDR-H2', 'CDR-H1', 'CDR-L3', 'CDR-L2', 'CDR-L1']
    for segment in segments:
        get_pattern_residue_decomposition('abdb_outfiles/respairs_epitope_segment_notationx.csv', 'epitope',
                                          'agresnumiset', segment)

def prepdata_pattern_interaction_map(abfile, agfile):
    '''
    Preps neat data for intearaction map
    :return:
    '''
    agfile = 'abdb_outfiles/respairs_epitope_segment_notationx.csv'
    abfile = 'abdb_outfiles/respairs_paratope_segment_notationx.csv'
    abdf, agdf = pd.read_csv(abfile), pd.read_csv(agfile)
    abdict = {}
    agdict = {}
    for i,row in abdf.iterrows():
        pdbsegment = row.pdbid + row.segment
        abdict[pdbsegment] = row
    for i2,row2 in agdf.iterrows():
        pdbsegment = row2.pdbid + row2.segment
        agdict[pdbsegment] = row2
    # print(len(abdict.keys()))V
    # print(len(agdict.keys()))
    pdbsegments = [item for item in abdict.keys() if item in agdict]
    abag_data = []
    for pdbsegment in pdbsegments:
        abrow  = abdict[pdbsegment]
        agrow = agdict[pdbsegment]
        pdbid = abrow.pdbid
        abchain = abrow.abchain
        segment = abrow.segment
        paratope = abrow.paratope
        plen = abrow.plen
        abresnumiset = abrow.abresnumiset
        ab_gap_pattern = abrow.gap_patterns
        absegment = abrow.absegment
        ab_gapstatus = abrow.gapstatus
        epitope = agrow.epitope
        elen = agrow.elen
        agresnumiset = agrow.agresnumiset
        ag_gap_pattern = agrow.gap_patterns
        ag_gapstatus = agrow.gapstatus
        agchain = agrow.agchain
        abag_datum = [pdbid, abchain, segment, paratope, plen, abresnumiset, ab_gap_pattern, absegment, ab_gapstatus,
                     epitope, elen, agresnumiset, ag_gap_pattern, ag_gapstatus, agchain]
        abag_data.append(abag_datum)
    columns = ['pdbid', 'abchain', 'segment', 'paratope', 'plen', 'abresnumiset', 'ab_gap_pattern', 'absegment',
               'ab_gapstatus', 'epitope', 'elen', 'agresnumiset', 'ag_gap_pattern', 'ag_gapstatus', 'agchain']
    outdf = pd.DataFrame(abag_data, columns=columns)
    print(outdf.head())
    outname = 'abdb_outfiles/respairs_segment_notationx_merged.csv'
    outdf.to_csv(outname, index=False)
    sys.exit()
    for i,row in outdf.iterrows():
        pdbdf = outdf[(outdf.pdbid==row.pdbid)]
        chains = pdbdf.abchain.unique()
        for chain in chains:
            chaindf = pdbdf[pdbdf.abchain == chain]
            paratope_cdr = ''
            paratope_fr = ''
            cdrs  = ['CDR-H1', 'CDR-H2', 'CDR-H3']
            for cdr in cdrs:
                cdrdf = chaindf[chaindf.segment==cdr]
                try:
                    paratope_cdr += 'J'# + cdrdf.ab_gap_pattern
                except:
                    paratope_cdr += 'j'
            print(paratope_cdr)
            sys.exit()

        sys.exit()
    sys.exit()
    outname = 'abdb_outfiles/respairs_segment_notationx_meVrged.csv'
    outdf.to_csv(outname, index=False)

def plot_pattern_interaction_map(infile, segment):
    '''
    plots gap pattern interaction map
    :param infile:
    :return:
    '''
    df = pd.read_csv(infile)
    # df = df[df.segment == 'CDR-H2']
    if segment != 'all':
        df = df[df.segment==segment]

    ctdf = pd.crosstab(df.ab_gap_pattern, df.ag_gap_pattern)
    ctdict = {'log_frequency':[], 'ab_gap_pattern': [], 'ag_gap_pattern': [], 'ab_ag_gap_pattern':[], 'frequency':[]}
    for i, row in ctdf.iterrows():
        frequencies = row.tolist()
        n = len(frequencies)
        ab_gap_pattern = [row.name]*n
        ag_gap_pattern = row.index.tolist()
        ab_ag_gap_pattern = ['%s-%s' % (ab,ag) for ab, ag in zip(ab_gap_pattern,ag_gap_pattern)]
        # print(len(frequencies), len(ab_gap_pattern), len(ag_gap_pattern))
        # ctdict['frequency'] += frequencies
        #ignores 0 frequency pairs
        for f,agp in zip(frequencies, ag_gap_pattern):
            if f > 0:
                ctdict['log_frequency'] += list(np.round(np.log10([f]),2))*f
                ctdict['frequency'] += [f]*f
                ctdict['ab_gap_pattern'] += [i]*f
                ctdict['ag_gap_pattern'] += [agp]*f
                ctdict['ab_ag_gap_pattern'] += ['%s-%s'% (i,agp)]*f
    topn = 10
    freqdf = pd.DataFrame(ctdict)
    heatdf = freqdf.drop_duplicates()
    topheatdf = heatdf.sort_values(by=['log_frequency']).iloc[:] #[-topn:]
    ab_ag_gap_pattern_list = freqdf.ab_ag_gap_pattern.value_counts().index.tolist()[:] # [:topn]
    topnfreqdf = freqdf[freqdf.ab_ag_gap_pattern.isin(ab_ag_gap_pattern_list)]
    ab_ag_gap_pattern_ord = CategoricalDtype(categories=ab_ag_gap_pattern_list, ordered=True)
    topnfreqdf['ab_ag_gap_pattern_ord'] = topnfreqdf['ab_ag_gap_pattern'].astype(str).astype(ab_ag_gap_pattern_ord)
    g = (ggplot(topheatdf)
         + aes('ab_gap_pattern', 'ag_gap_pattern', fill='log_frequency')
         + geom_tile()
         + geom_text(aes(label='log_frequency'), color='white')
         + coord_equal()
         + theme(axis_text_x=element_text(rotation=90))
        )
    g2 = (ggplot(topnfreqdf)
          + aes(x = 'ab_ag_gap_pattern_ord')
          + geom_bar()
          # + theme(axis_text_x=element_text(rotation=90))
          + coord_flip()
         )
    # gview_temp(g)
    # gview_temp(g2)
    outfigbar = infile.split('/')[-1].split('.')[0] + '_top%s_pattern_pairs_bar_%s' % (topn, segment)
    outfigheatmap = infile.split('/')[-1].split('.')[0] + '_top%s_pattern_pairs_heatmap_%s' % (topn, segment)
    # gsave_pdf(g2,outfigbar)
    # gsave_pdf(g,outfigheatmap)
    return topnfreqdf, topn


def plot_pattern_interaction_map_all_segments():
    '''
    plots pattern interaction map for all segments
    outputs *_intermap.csv
    :return:
    '''
    segments = ['CDR-H3', 'CDR-H2', 'CDR-H1', 'CDR-L3', 'CDR-L2', 'CDR-L1', 'all', 'HFR1', 'HFR2', 'HFR3', 'LFR1',
                'LFR2', 'LFR3', 'HFR4']
    segmentdfs = []
    for segment in segments:
        topnfreqdf, topn = plot_pattern_interaction_map('abdb_outfiles/respairs_segment_notationx_merged.csv', segment)
        nrows = topnfreqdf.shape[0]
        topnfreqdf['segment'] = [segment]*nrows
        absegment = re.sub('[LH-]','',segment)
        if segment.startswith('C'):
            abchain = segment[-2]
        elif segment[1] == 'F':
            abchain = segment[0]
        topnfreqdf['abchain']  = [abchain]*nrows
        topnfreqdf['absegment'] = [absegment]*nrows
        sorted_topnfreqdf = topnfreqdf.sort_values('frequency', ascending=False).drop_duplicates(
            subset='ab_ag_gap_pattern')
        ranks = [i for i in range(1,sorted_topnfreqdf.shape[0]+1)]
        sorted_topnfreqdf['order'] = ranks
        segmentdfs.append(sorted_topnfreqdf)
    intermapdf = pd.concat(segmentdfs)
    print(intermapdf.head())
    # intermapdf_unique = intermapdf.drop_duplicates('ab_ag_gap_pattern')
    # g = (ggplot(intermapdf_unique)
    g = (ggplot(intermapdf)
         + aes('order', 'frequency')
         # + aes('ab_ag_gap_pattern_ord')
         # + geom_col()
         # + geom_bar(stat='identity')
         + geom_jitter()
         + coord_flip()
         # + facet_grid('absegment~abchain', scales='free')
         + facet_wrap('segment', scales='free', ncol=3)
         + geom_text(aes(label='ab_ag_gap_pattern'),  # new
                     size=2, va='bottom', nudge_x=0, nudge_y=5 )
         + theme(figure_size=(25,25), panel_spacing=0.5)
        )
    gview_temp(g)
    sys.exit()
    outfig = 'intermap_top%s_pattern' % topn
    gsave_pdf(g,outfig)

def plot_residue_interaction_map(infile, segment):
    '''
    preps neat data for resisdue map interaction
    :param infille:
    :return:
    '''
    df = pd.read_csv(infile)
    df = df[df.segment == segment]
    ctdf = pd.crosstab(df.abres, df.agres)
    ctdict = {'frequency':[], 'abres':[], 'agres':[]}
    ncol = ctdf.shape[-1]
    print(df.head())
    for i,row in ctdf.iterrows():
        abres = [row.name]*ncol
        agres = row.index.tolist()
        frequencies = list(row.values)
        ctdict['frequency'] += frequencies
        ctdict['abres'] += abres
        ctdict['agres'] += agres
    freqdf = pd.DataFrame(ctdict)
    g = (ggplot(freqdf)
         + aes('abres', 'agres', fill='frequency')
         + geom_tile()
         + coord_equal()
         + geom_text(aes(label='frequency'), color='white')
         + theme(figure_size=(10,10))
    )
    # gview_temp(g)
    # outfig = infile.split('/')[-1].split('.')[0] + '_intermap_residue_%s'%segment
    # gsave_pdf(g,outfig)
    # visualize as a network
    plt.figure(figsize=(25,25))
    G = nx.Graph()
    G.clear()
    node_colors = []
    for i,row in freqdf.iterrows():
        agres = row.agres + '+'
        if row.frequency > 0:
            G.add_edge(row.abres, agres, weight=row.frequency)
    for node in G:
        if '+' in node:
            node_colors.append('orange')
        else:
            node_colors.append('salmon')
    edges = G.edges()
    weights = [G[u][v]['weight']**0.8 for u,v in edges]
    node_sizes = [5000]*len(node_colors)
    nx.draw_spring(G,edges=edges, width=weights, with_labels=True, alpha=0.9, node_color=node_colors,
                 node_size=node_sizes,
            font_color='white', font_weight='bold', font_size=20)
    outfig2 = infile.split('/')[-1].split('.')[0] + '_intermap_residue_network_%s'%segment
    metrics = []
    # print(sorted(nx.degree(G), key= lambda item: item[1], reverse=True)[:3])
    degree = [('degree')] + ['%s: %s' % (res,round(val, 2)) for res, val in sorted(nx.degree_centrality(
        G).items(), key=lambda item: item[1],  reverse=True)[:]]
    betweenness = [('betweenness')] + ['%s: %s' % (res,round(val, 2)) for res, val in sorted(nx.betweenness_centrality(
        G).items(), key=lambda item: item[1],  reverse=True)[:]]
    closeness = [('closeness')] + ['%s: %s' % (res,round(val, 2)) for res, val in sorted(nx.closeness_centrality(
        G).items(), key=lambda item: item[1],  reverse=True)[:]]
    n_nodes = len(node_colors)
    density = ['density'] + [round(nx.density(G),3)]*n_nodes
    for item in [degree, betweenness,closeness,density]:
        filetype = infile.split('_')[2]
        metric = [filetype, segment] + item
        print(filetype)
        metrics.append(metric)
    metric_columns = ['origin', 'segment', 'metric']
    metric_columns += ['residue%s'% i for i in range(1,n_nodes+1)]
    print(metric_columns)
    metricdf  = pd.DataFrame(metrics,columns= metric_columns)
    metrics2 = []
    for zitems in zip(metric_columns[3:], degree[3:], betweenness[3:],closeness[3:]):
        metrics2.append(zitems)
        print(zitems)
    metric_columns2 = ['residue', 'degree', 'betweenness', 'closeness']
    metricdf2 = pd.DataFrame(metrics2, columns = metric_columns2)
    outfile = 'abdb_outfiles/' + infile.split('/')[-1].split('.')[0] + '_intermap_residue_network_stats_%s.csv'%segment
    metricdf2.to_csv(outfile, index=False)
    # sys.exit()
    save_pdf(outfig2)

#
def plot_residue_position_interaction_map(infile, segment, abchain):
    '''
    plots residue-position interactin map.
    :return:
    '''
    df = pd.read_csv(infile)
    # df = df[df.abchain == abchain]
    if segment !='all':
        df = df[df.segment == segment]
    if abchain != 'all':
        df = df[df.abchain == abchain]
    # abrespos = ['%s%s' %(row.abres, row.abresnum) for i, row in df.iterrows()]
    abrespos = ['%s%s' %(row.abresnum, row.abres) for i, row in df.iterrows()]
    abpos = ['%s' % row.abresnum for i, row in df.iterrows()]
    agrespos = ['%s%s' %(row.agres, row.agresnum) for i, row in df.iterrows()]
    df['abrespos'] = abrespos
    df['agrespos'] = agrespos
    df['segment_cat'] = df['segment'].astype('category')
    ctdf = pd.crosstab(df.abrespos, df.agres)
    ctdict = {'frequency':[], 'abrespos':[], 'agres':[], 'abpos':[]}
    ncol = ctdf.shape[-1]
    for i,row in ctdf.iterrows():
        abrespos = [row.name]*ncol
        abpos = [row.name[:-3]]*ncol
        agres = row.index.tolist()
        frequencies = list(row.values)
        ctdict['frequency'] += frequencies
        ctdict['abrespos'] += abrespos
        ctdict['agres'] += agres
        ctdict['abpos'] += abpos
    freqdf = pd.DataFrame(ctdict)
    percent = 1
    n = int(freqdf.shape[0]*percent) # get percent from total
    percentstr = 'top%spercent' % str(int(percent*100))
    topnrespos = freqdf.sort_values(by='frequency').iloc[-n:]
    # g = (ggplot(topnrespos)
    #      + aes('agres', 'abrespos', fill='frequency')
    #      + geom_tile()
    #      + coord_equal()
    #      + theme(figure_size=(10,10))
    #      + geom_text(aes(label='frequency'), color='white')
    # )
    # outfig = infile.split('/')[-1].split('.')[0] + '_intermap_residue_position_abchain%s_segment%s_%s' % (
    #     abchain, segment, percentstr)
    # gsave_pdf(g,outfig)
    topnrespos_dict_ag  = {}
    for i, row in topnrespos.iterrows():
        if row.abpos not in topnrespos_dict_ag:
            topnrespos_dict_ag[row.abpos] = {row.agres:row.frequency}
        elif row.agres not in topnrespos_dict_ag[row.abpos]:
            topnrespos_dict_ag[row.abpos][row.agres] = row.frequency
        else:
            topnrespos_dict_ag[row.abpos][row.agres] += row.frequency
    topnrespos_dict_ab = {}
    for iab, rowab in topnrespos.iterrows():
        abres = rowab.abrespos[-3:]
        if rowab.abpos not in topnrespos_dict_ab:
            topnrespos_dict_ab[rowab.abpos]={abres: rowab.frequency}
        elif abres not in topnrespos_dict_ab[rowab.abpos]:
            topnrespos_dict_ab[rowab.abpos][abres] = rowab.frequency
        else:
            topnrespos_dict_ab[rowab.abpos][abres] += rowab.frequency
    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    aas = sorted(aadict.items(), key= lambda item: item[-1])
    topnrespos_pwm_ag = []
    positions_ag = natsort.natsorted(topnrespos_dict_ag)
    for position in positions_ag:
        datum = [position]
        for aa in aas:
            try:
                aa_freq = topnrespos_dict_ag[position][aa[0]]
                datum.append(aa_freq)
            except KeyError:
                datum.append(0)
        topnrespos_pwm_ag.append(datum)
    columns = ['P0'] + [item[1] for item in aas]
    # pwmdf = pd.DataFrame(topnrespos_pwm_ag,columns=columns)
    outlist = [columns] + topnrespos_pwm_ag
    outstr = ''
    for line in outlist:
        str_line = '\t'.join([str(item) for item in line]) + '\n'
        outstr += str_line
    outname_ag = 'abdb_outfiles/' + infile.split('/')[-1].split('.')[0] + '_pwm_ag_abchain%s_segment%s_%s.tsv' % (
               abchain, segment, percentstr)
    outfile = open(outname_ag, 'w')
    outfile.write(outstr)
    outfile.close()
    pos_annotations = ','.join(positions_ag)
    weblogo_outpdf = 'abdb_figures/' + outname_ag.split('/')[-1].split('.')[0] + '_weblogo.pdf'
    weblogo_args = (outname_ag, 'transfac', weblogo_outpdf, 'pdf', 'protein', 'yes', pos_annotations,'probability',
                    'yes', segment, 'chemistry')
    weblogo_command = 'weblogo -f %s -D %s -o %s -F %s -A %s --scale-width %s --annotate %s -U %s --rotate-number ' \
                      '%s -P %s -c %s' % \
                      weblogo_args
    print(weblogo_command)
    os.popen(weblogo_command)
    os.system('open %s' % weblogo_outpdf)
    # start pwm ab bits
    topnrespos_pwm_ab = []
    positions_ab = natsort.natsorted(topnrespos_dict_ab)
    for position in positions_ab:
        datum = [position]
        for aa in aas:
            try:
                aa_freq = topnrespos_dict_ab[position][aa[0]]
                datum.append(aa_freq)
            except KeyError:
                datum.append(0)
        topnrespos_pwm_ab.append(datum)
    outlist_ab = [columns] + topnrespos_pwm_ab
    outstr_ab = ''
    for line in outlist_ab:
        str_line = '\t'.join([str(item) for item in line]) + '\n'
        outstr_ab += str_line
    outname_ab = 'abdb_outfiles/' + infile.split('/')[-1].split('.')[0] + '_pwm_ab_abchain%s_segment%s_%s.tsv' % (
        abchain, segment, percentstr)
    outfile_ab = open(outname_ab, 'w')
    outfile_ab.write(outstr_ab)
    outfile_ab.close()
    pos_annotations = ','.join(positions_ab)
    weblogo_outpdf_ab = 'abdb_figures/' + outname_ab.split('/')[-1].split('.')[0] + '_weblogo.pdf'
    weblogo_args_ab = (outname_ab, 'transfac', weblogo_outpdf_ab, 'pdf', 'protein', 'yes', pos_annotations,
                       'probability',
                    'yes', segment, 'chemistry')
    weblogo_command_ab = 'weblogo -f %s -D %s -o %s -F %s -A %s --scale-width %s --annotate %s -U %s --rotate-number ' \
                      '%s -P %s -c %s' % \
                      weblogo_args_ab
    os.popen(weblogo_command_ab)
    os.system('open %s' % weblogo_outpdf_ab)
    print(outstr)
    print(outstr_ab)

def plot_gap_pattern_pwm(infile, pattern, segment):
    '''
    plots residue 'pairwise' pwm to see sort of long range relationship amongs residues
    :return:
    '''
    df = pd.read_csv(infile)
    df = df[(df.gap_patterns == pattern) & (df.segment== segment)]
    # make fasta file
    content = ''
    for i,row in df.iterrows():
        try:
            residues = row.paratope
        except AttributeError:
            residues = row.epitope
        for i2 in range(len(residues)-1):
            # seq = '>\n' + '-'*i2 + residues[i2]
            seq = '-'*i2 + residues[i2]
            for i3, residue in enumerate(residues[i2+1:]):
                seq2 = seq + '-'*i3 + residue
                seq3 =  '>\n' + seq2 + '-'*(len(residues)-len(seq2)) + '\n'
                content += seq3
    outname = 'abdb_outfiles/' + infile.split('/')[-1].split('.')[0] + '_%s_%s.fasta' % (pattern, segment)
    outfile = open(outname, 'w')
    outfile.write(content)
    outfile.close()
    weblogo_outpdf = 'abdb_figures/' + outname.split('/')[-1].split('.')[0] + '_weblogo.pdf'
    print(weblogo_outpdf)
    # pos_annotations = ','.join([str(item) for item in range(1,len(pattern)+1)])
    pos_annotations = ','.join([item for item in pattern])
    weblogo_args = (outname, 'fasta', weblogo_outpdf, 'pdf', 'protein', 'yes', pos_annotations,
                       'probability', 'no', '_', 'chemistry', '6', '4')
    weblogo_command = 'weblogo -f %s -D %s -o %s -F %s -A %s --scale-width %s --annotate %s -U %s --rotate-number ' \
                         '%s -P %s -c %s --fontsize %s --number-fontsize %s' % \
                         weblogo_args
    os.system(weblogo_command)
    os.system('open %s' % weblogo_outpdf)


def loop_plot_gap_pattern_pwm():
    '''
    loops plot_gap_pattern_pwm
    :return:
    '''
    patterns =['XX', 'XXX', 'XXXX', 'XXXXX', 'XXXXXX', 'XXXXXXX', 'XXXXXXXXX']
    for pattern in patterns:
        plot_gap_pattern_pwm('abdb_outfiles/respairs_paratope_segment_notationx.csv', pattern, 'CDR-H3')

def get_context_overview():
    '''
    gets context from victor's drawing.
    :return:
    '''
    infile = 'abdb_outfiles/respairs_paratope_segment_notationx.csv'
    df = pd.read_csv(infile)
    # dfx = df[(df.gap_patterns == 'XX') & (df.segment == 'CDR-H3')]
    # dfx = df[(df.gap_patterns == 'XXXXXXX')]# & (df.segment == 'CDR-H3')]
    dfx = df.dropna()
    edges = []
    edges_next = []
    edges_next_single = []
    for i,row in dfx.iterrows():
        paratope = row.paratope
        lenp = len(paratope)
        for i2, res in enumerate(paratope):
            next_i = -(lenp-i2-1)
            next_res = paratope[next_i:]
            prev_res = paratope[:i2]
            if next_res == paratope:
                context = prev_res
            else:
                context = prev_res + next_res
                context_next = next_res
                context_next_single = next_res[0]
            edge = (res, context)
            edge_next = (res, context_next)
            edge_next_single = (res,context_next_single)
            edges.append(edge)
            edges_next.append(edge_next)
            edges_next_single.append(edge_next_single)
    count_dict = {}
    for edge in edges_next_single:
        if edge not in count_dict:
            count_dict[edge] = 1
        else:
            count_dict[edge] += 1
    G = nx.Graph()
    plt.figure(figsize=(15,15))
    edges_weights = count_dict.items()
    # edges_weights = [item for item in edges_weights if item[0][0] == 'D']
    # print(edges_weights)
    for edge, weight in edges_weights:
        G.add_edge(edge[0],edge[1],weight=weight)
    edges = G.edges()
    weights = [G[u][v]['weight'] for u,v in edges]
    nx.draw(G,edges=edges,width=weights, with_labels=True, alpha=0.3, font_color='black')
    pdfname = 'abdb_figures/temp_net.pdf'
    outname = 'context_single_overview'
    # plt.savefig(pdfname)
    # os.system('open %s' %pdfname)
    save_pdf(outname)

def get_context_pattern(infile,segment):
    '''
    gets context from victor's drawing.
    :return:
    '''
    # infile = 'abdb_outfiles/respairs_paratope_segment_notationx.csv'
    df = pd.read_csv(infile)
    dfx = df.dropna()
    # patterns = ['XX', 'XXX', 'XXXX', 'XXXXX']
    patterns = ['XXXXX']
    if 'epitope' in infile:
        patterns = ['XX', 'XXXXX']
    nrow,ncol  = 2,2
    for pattern_i, pattern in enumerate(patterns):
        plt.figure(figsize=(15,15))
        coord = int('%s%s%s' % (nrow,ncol,pattern_i+1))
        # plt.subplot(coord)
        # ax = plt.gca()
        # ax.set_xlabel(pattern)
        # ax.set_yticklabels([])
        # ax.set_xticklabels([])
        # ax.tick_params(axis=u'both', which=u'both', length=0)
        # # ax.yaxis.set_visible(False)
        # # ax.xaxis.set_visible(False)
        # # ax.set_xlabel(pattern)
        # dfx = df[(df.gap_patterns == pattern)]# & (df.segment == 'CDR-H3')]
        dfx = df[(df.gap_patterns == pattern) & (df.segment == segment)]
        edges = []
        edges_next = []
        edges_next_single = []
        edges_prev_next_single = []
        for i,row in dfx.iterrows():
            try: # handles paratope and epitope
                paratope = row.paratope
            except AttributeError:
                paratope = row.epitope
            lenp = len(paratope)
            for i2, res in enumerate(paratope):
                next_i = -(lenp-i2-1)
                next_res = paratope[next_i:]
                prev_res = paratope[:i2]
                if i2 == 0:
                    prev_res = '-'
                if next_res == paratope: # check for last residue
                    context = prev_res
                    context_next_single = '-'
                    context_prev_next_single = prev_res[-1] + '-'
                else:
                    context = prev_res + next_res
                    context_next = next_res
                    context_next_single = next_res[0]
                    context_prev_next_single = prev_res[-1] + next_res[0]
                    if res == context_next_single:
                        context_next_single += "'"
                edge = (res, context)
                edge_next = (res, context_next)
                edge_next_single = (res, context_next_single)
                edge_prev_next_single = (res, context_prev_next_single)
                edges.append(edge)
                edges_next.append(edge_next)
                if i2 != lenp-1:
                    edges_next_single.append(edge_next_single)
                edges_prev_next_single.append(edge_prev_next_single)
                # print(paratope, res, context_prev_next_single, context_next_single)
        tag = infile.split('/')[-1].split('.')[0]
        outnames = [tag + '_context_next_single', tag + '_context_prev_next_single']
        context_edges = [edges_next_single, edges_prev_next_single]
        for outname, context_edge in zip(outnames[:1],context_edges[:1]):
            count_dict = {}
            for edge in context_edge:
                if edge not in count_dict:
                    count_dict[edge] = 1
                else:
                    count_dict[edge] += 1
            G = nx.Graph()
            edges_weights = sorted(count_dict.items(), key= lambda item: item[-1])
            for edge, weight in edges_weights:
                G.add_edge(edge[0],edge[1],weight=weight)
            edges = G.edges()
            weights = [G[u][v]['weight'] for u,v in edges]
            node_size = [5000]*len(weights)
            font_size = 50
            metrics = []
            # print(sorted(nx.degree(G), key= lambda item: item[1], reverse=True)[:3])
            degree = [('degree')] + ['%s: %s' % (res,round(val, 2)) for res, val in sorted(nx.degree_centrality(
                G).items(), key=lambda item: item[1],  reverse=True)[:3]]
            betweenness = [('betweenness')] + ['%s: %s' % (res,round(val, 2)) for res, val in sorted(nx.betweenness_centrality(
                G).items(), key=lambda item: item[1],  reverse=True)[:3]]
            closeness = [('closeness')] + ['%s: %s' % (res,round(val, 2)) for res, val in sorted(nx.closeness_centrality(
                G).items(), key=lambda item: item[1],  reverse=True)[:3]]
            density = ['density'] + [round(nx.density(G),3)]*3
            for item in [degree, betweenness,closeness,density]:
                filetype = infile.split('_')[2]
                metric = [filetype, segment, pattern] + item
                print(filetype)
                metrics.append(metric)
            metric_columns = ['origin', 'segment', 'pattern', 'metric', 'residue1: score', 'residue2: score',
                            'residue3: score']
            metricdf  = pd.DataFrame(metrics,columns= metric_columns)
            outmetric_name = 'abdb_outfiles/'+ outname + '_network_metrics_%s.csv' % pattern
            metricdf.to_csv(outmetric_name, index=False)
            # print(metricdf)
            # sys.exit()
            if 'prev_next' in outname:
                node_size = [200]* len(weights)
                font_size = 10
            nx.draw_spring(G,edges=edges,width=weights, with_labels=True, alpha=0.5, font_color='black',
            node_size=node_size,
                    font_size=font_size, font_weight='bold')
            outname += '_%s_%s' % (pattern,segment)
            save_pdf(outname)
            outgraph = 'abdb_outfiles/' + outname.split('.')[0] + '.ssv'
            print(outgraph)
            # nx.write_graphml(G,outgraph)
            nx.write_weighted_edgelist(G,outgraph)
            G.clear()
            # sys.exit()


def get_context_pattern_top3(infile):
    '''
    gets context from top 3 epitope motifs.
    for single motif manually for now
    :return:
    '''
    # infile = 'abdb_outfiles/respairs_paratope_segment_notationx.csv'
    df = pd.read_csv(infile)
    dfx = df.dropna()
    top3 = dfx.gap_patterns.value_counts()[:3]
    patterns = top3.index
    # patterns = ['XXX', 'XXXX', 'XXXXX'] #paratope
    # patterns = ['XXX']
    # patterns = ['XXXX']
    # patterns = ['XXXXX']
    ptag = '_'.join(patterns)
    dfxs = []
    for pattern in patterns:
        dfp = dfx[df.gap_patterns == pattern]
        dfxs.append(dfp)
    dfxc = pd.concat(dfxs)
    edges = []
    edges_next = []
    edges_next_single = []
    edges_prev_next_single = []
    for i,row in dfxc.iterrows():
        if 'epitope' in infile:
            paratope = row.epitope
        elif 'paratope' in infile:
            paratope = row.paratope
        lenp = len(paratope)
        for i2, res in enumerate(paratope):
            next_i = -(lenp-i2-1)
            next_res = paratope[next_i:]
            prev_res = paratope[:i2]
            if i2 == 0:
                prev_res = '-'
            if next_res == paratope: # check for last residue
                context = prev_res
                context_next_single = '-'
                context_prev_next_single = prev_res[-1] + '-'
            else:
                context = prev_res + next_res
                context_next = next_res
                context_next_single = next_res[0]
                context_prev_next_single = prev_res[-1] + next_res[0]
                if res == context_next_single:
                    context_next_single += "'"
            edge = (res, context)
            edge_next = (res, context_next)
            edge_next_single = (res, context_next_single)
            edge_prev_next_single = (res, context_prev_next_single)
            edges.append(edge)
            edges_next.append(edge_next)
            if i2 != lenp-1:
                edges_next_single.append(edge_next_single)
            edges_prev_next_single.append(edge_prev_next_single)
            print(paratope, res, context_prev_next_single, context_next_single)
        print('%s/%s: %s' % (i, dfxc.shape[0],infile))
    edge_next_count = {}
    for pair in edges_next_single:
        if pair not in edge_next_count:
            edge_next_count[pair] = 1
        else:
            edge_next_count[pair] += 1
    edge_nextdf = pd.DataFrame([[item[0][0],item[0][1], item[1]] for item in edge_next_count.items()],
                               columns=['source', 'target', 'count']).sort_values(by='count')
    if 'paratope' in infile:
        tag = 'paratope'
    elif 'epitope' in infile:
        tag = 'epitope'
    if 'dewitt' in infile:
        tag = infile.split('/')[-1].split('.')[0]
    outnext = 'abdb_outfiles/%s_%s_top3motif_respairs_next.csv' % (tag,ptag)
    edge_nextdf.to_csv(outnext, index=False)
    # sys.exit()
    # edge_prev_next_count = {}
    # for pair in edges_prev_next_single:
    #     if pair not in edge_prev_next_count:
    #         edge_prev_next_count[pair] = 1
    #     else:
    #         edge_prev_next_count[pair] += 1
    # edge_prev_nextdf = pd.DataFrame([[item[0][0],item[0][1], item[1]] for item in edge_prev_next_count.items()],
    #                            columns=['source', 'target', 'count']).sort_values(by='count')
    # outprevnext = 'abdb_outfiles/top3motif_respairs_prev_next.csv'
    # edge_prev_nextdf.to_csv(outprevnext, index=False)
    # print(sum(edge_nextdf['count']))
    # print(sum(edge_prev_nextdf['count']))
    # sys.exit()


def get_context_pattern_top3_position(infile):
    '''
    gets context from top 3 epitope motifs.
    for single motif manually for now
    account for position
    :return:
    '''
    # infile = 'abdb_outfiles/respairs_paratope_segment_notationx.csv'
    df = pd.read_csv(infile)
    dfx = df.dropna()
    top3 = dfx.gap_patterns.value_counts()[:3]
    patterns = top3.index
    # patterns = ['XXX', 'XXXX', 'XXXXX'] #paratope
    # patterns = ['XXX']
    # patterns = ['XXXX']
    # patterns = ['XXXXX']
    ptag = '_'.join(patterns)
    dfxs = []
    for pattern in patterns:
        dfp = dfx[df.gap_patterns == pattern]
        dfxs.append(dfp)
    dfxc = pd.concat(dfxs)
    edges = []
    edges_next = []
    edges_next_single = []
    edges_prev_next_single = []
    for i,row in dfxc.iterrows():
        if 'epitope' in infile:
            paratope = row.epitope
        elif 'paratope' in infile:
            paratope = row.paratope
        lenp = len(paratope)
        for i2, res in enumerate(paratope):
            next_i = -(lenp-i2-1)
            next_res = paratope[next_i:]
            prev_res = paratope[:i2]
            if i2 == 0:
                prev_res = '-'
            if next_res == paratope: # check for last residue
                context = prev_res
                context_next_single = '-'
                context_prev_next_single = prev_res[-1] + '-'
            else:
                context = prev_res + next_res
                context_next = next_res
                context_next_single = next_res[0]
                context_prev_next_single = prev_res[-1] + next_res[0]
                if res == context_next_single:
                    context_next_single += "'"
            edge = (res, context)
            edge_next = (res, context_next)
            # edge_next_single = (res, context_next_single)
            # add resindex as position
            edge_next_single = (res+str(row.resindex), context_next_single+str(row.resindex))
            edge_prev_next_single = (res, context_prev_next_single)
            edges.append(edge)
            edges_next.append(edge_next)
            if i2 != lenp-1:
                edges_next_single.append(edge_next_single)
            edges_prev_next_single.append(edge_prev_next_single)
            print(paratope, res, context_prev_next_single, context_next_single)
        print('%s/%s: %s' % (i, dfxc.shape[0],infile))
    edge_next_count = {}
    for pair in edges_next_single:
        if pair not in edge_next_count:
            edge_next_count[pair] = 1
        else:
            edge_next_count[pair] += 1
    edge_nextdf = pd.DataFrame([[item[0][0],item[0][1], item[1]] for item in edge_next_count.items()],
                               columns=['source', 'target', 'count']).sort_values(by='count')
    if 'paratope' in infile:
        tag = 'paratope'
    elif 'epitope' in infile:
        tag = 'epitope'
    if 'dewitt' in infile:
        tag = infile.split('/')[-1].split('.')[0]
    outnext = 'abdb_outfiles/%s_%s_top3motif_respairs_next_position.csv' % (tag,ptag)
    edge_nextdf.to_csv(outnext, index=False)
    # sys.exit()
    # edge_prev_next_count = {}
    # for pair in edges_prev_next_single:
    #     if pair not in edge_prev_next_count:
    #         edge_prev_next_count[pair] = 1
    #     else:
    #         edge_prev_next_count[pair] += 1
    # edge_prev_nextdf = pd.DataFrame([[item[0][0],item[0][1], item[1]] for item in edge_prev_next_count.items()],
    #                            columns=['source', 'target', 'count']).sort_values(by='count')
    # outprevnext = 'abdb_outfiles/top3motif_respairs_prev_next.csv'
    # edge_prev_nextdf.to_csv(outprevnext, index=False)
    # print(sum(edge_nextdf['count']))
    # print(sum(edge_prev_nextdf['count']))
    # sys.exit()


def merge_segment_motif(infile):
    '''
    merge motif in each segment
    :return:
    '''
    df = pd.read_csv(infile)
    # print(df.head())
    # sys.exit()
    chains = df.abchain.unique()
    hsegments= ['HFR1', 'CDR-H1', 'HFR2', 'CDR-H2', 'HFR3', 'CDR-H3', 'HFR4']
    lsegments= ['LFR1', 'CDR-L1', 'LFR2', 'CDR-L2', 'LFR3', 'CDR-L3', 'LFR4']
    merged_motifs = []
    merged_chains = []
    merged_pdbids = []
    merged_pmotifs = []
    unique_pdbids =df.pdbid.unique()
    for pdbid in unique_pdbids:
        pdbdf = df[df.pdbid == pdbid]
        for chain in chains:
            if chain == 'L':
                segments = lsegments
            else:
                segments = hsegments
            merged_motif = ''
            merged_pmotif = ''
            for segment in segments:
                segdf = pdbdf[pdbdf.segment == segment]
                if segdf.shape[0] != 0:
                    smotif = segdf.gap_patterns.values[0] + 'j'
                    pmotif = segdf.abresnumiset.values[0] + 'j'
                else:
                    smotif = 'j'
                    pmotif = 'j'
                merged_motif += smotif
                merged_pmotif += pmotif
            merged_motifs.append(merged_motif)
            merged_pmotifs.append(merged_pmotif)
            merged_chains.append(chain)
            merged_pdbids.append(pdbid)
    merged_df = pd.DataFrame()
    merged_df['pdbid'] = merged_pdbids
    merged_df['merged_motif'] = merged_motifs
    merged_df['merged_pmotif'] = merged_pmotifs
    merged_df['abchain'] = merged_chains
    outfile = infile.split('.')[0] + '_motif_merged.csv'
    merged_df.to_csv(outfile, index=False)

def get_levenshtein_distance(infile):
    '''
    gets pairwise levenshtein distance from motifs
    :param infile:
    :return:
    '''
    df = pd.read_csv(infile)
    df = df[df.abchain=='H']
    distances = []
    for mm in df.merged_motif:
        distance = []
        for mm2 in df.merged_motif:
            ld = jellyfish.levenshtein_distance(mm,mm2)
            distance.append(ld)
        distances.append(distance)
    print(distances)
    ddf = pd.DataFrame(distances)
    print(df)

def add_motif_stats(infile):
    '''
    get motif limits
    :param infile:
    :return:
    '''
    df = pd.read_csv(infile)
    freq_allsegments = []
    freq_segments = []
    motif_lengths = []
    # df = df[(df.gap_patterns == 'XXXXX')]
    # print(df)
    # sys.exit()
    for i, row in df.iterrows():
        motif = row.gap_patterns
        mdf = df[df.gap_patterns == motif]
        freq_allsegment = mdf.shape[0]
        segment = row.segment
        segdf = mdf[mdf.segment == segment]
        freq_segment = segdf.shape[0]
        print(freq_allsegment,freq_segment)
        freq_segments.append(freq_segment)
        freq_allsegments.append(freq_allsegment)
        print(row)
        if 'epitope' in infile:
            if row.elen > 1:
                motif_len = row.elen + len([item for item in row.egapset.split('-') if item != '0'])
            else:
                motif_len = row.elen
        else:
            if row.plen > 1:
                motif_len = row.plen + len([item for item in row.gapset.split('-') if item != '0'])
            else:
                motif_len = row.plen
        motif_lengths.append(motif_len)
    df['motif_freq_allsegment'] = freq_allsegments
    df['motif_freq_segment'] = freq_segments
    df['motif_len'] = motif_lengths
    print(df.sort_values(by= 'motif_len', ascending=False).head())
    df.to_csv(infile, index=False)
    sys.exit()

def get_full_segment_seq(infile):
    '''
    gets full sequence instead of just intearcting residues.
    :param infile:
    :return:
    '''
    df = pd.read_csv(infile).iloc[:]
    data = []
    counter = 0
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
        try:
            amolecule = pdbcontent[-2][-1].strip()
            aspecies = pdbcontent[-1][-1].strip()
            lmolecule = pdbcontent[0][-1].strip()
            lspecies = pdbcontent[1][-1].strip()
            hmolecule = pdbcontent[2][-1].strip()
            hspecies = pdbcontent[3][-1].strip()
        except:
            amolecule = '-'
            aspecies = '-'
            lmolecule = '-'
            lspecies = '-'
            hmolecule = '-'
            hspecies = '-'
            counter += 1
            print(counter)
            print(open(pdbfile).read()[:600])
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
                datum = [pdbid, seq2, segment2, lvseq, lmolecule, lspecies,hmolecule,hspecies,amolecule,aspecies]
            elif 'H' in segment2:
                datum = [pdbid, seq2, segment2, hvseq, lmolecule, lspecies,hmolecule,hspecies,amolecule,aspecies]
            data.append(datum)
    columns = ['pdbid', 'segment_seq', 'segment', 'vgene', 'lmolecule', 'lspecies', 'hmolecule', 'hspecies',
               'amolecule', 'aspecies']
    df = pd.DataFrame(data, columns=columns)
    outfile = infile[:-4] + '_seq.csv'
    df.to_csv(outfile, index=False)
    print(outfile)
    print(df[df.segment=='CDR-H3'].segment_seq.value_counts())
    print(df[df.segment=='CDR-H3'].vgene.value_counts())

def get_vgenes():
    '''
    annotate vgenes using IMGT reference
    NOTE: using human and mouse only!
    :return:
    '''
    ighvfile = '../datasets/imgt_reference/IGHV.fasta'
    ighvfile = '../datasets/imgt_reference/IGHV.fasta'
    igvfiles = ['../datasets/imgt_reference/IGHV.fasta', '../datasets/imgt_reference/IGHV_mouse.fasta',
                '../datasets/imgt_reference/IGLV.fasta', '../datasets/imgt_reference/IGLV_mouse.fasta',
                '../datasets/imgt_reference/IGKV.fasta']
    # ighv = open(ighvfile).read().split('>')
    ighv = sum([open(igvfile).read().split('>') for igvfile in igvfiles], [])
    data = []
    specieses = ['Homo sapiens', 'Mus musculus_C57BL/6']
    functionalities = ['F', 'ORF', 'P']
    for v in ighv[:]:
        try:
            fastalines = v.splitlines()
            headers = fastalines[0].split('|')
            seqs = fastalines[1:]
            seq = ''.join(seqs).replace('.','')
            aaseq = Seq(seq).translate().tostring()
            functionality = headers[3]
            species = headers[2]
            imgt_genename = headers[1].split('*')[0]
            imgt_allelename = headers[1]
            chain = imgt_genename[2]
            ntcount = len(seq)
            triplet_status = ntcount%3
            if species in specieses and functionality in functionalities:
                datum = [seq, imgt_genename,imgt_allelename,functionality,species, chain, ntcount, triplet_status, aaseq]
                data.append(datum)
        except:
            print([v])
    columns = ['nt_seq', 'imgt_genename', 'imgt_allelename', 'functionality', 'species', 'chain', 'ntcount',
               'modulo_3', 'aa_seq']
    df = pd.DataFrame(data, columns=columns)
    df.to_csv('abdb_outfiles/imgt_vgenes.csv', index=False)
    print(df.head())
    sys.exit()

def annotate_vgenes():
    '''
    annoatate the pdb structures with the nearest vgenes (human and mouse only)
    :return:
    '''
    infile = 'abdb_outfiles/abdb_segment_seq.csv'
    vgene_file = 'abdb_outfiles/imgt_vgenes.csv'
    df = pd.read_csv(infile)
    vgenedf = pd.read_csv(vgene_file)
    available_species = ['HOMO SAPIENS', 'MUS MUSCULUS']
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
        if species not in available_species:
            species = 'MUS MUSCULUS'
        if species == 'HOMO SAPIENS':
            species = 'Homo sapiens' #IMGT formating
        elif species == 'MUS MUSCULUS':
            species = 'Mus musculus_C57BL/6'
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
    print(df.head())
    outfile = infile[:-4] + '_imgt_vgene.csv'
    print(outfile)
    df.to_csv(outfile, index=False)

def clean_annotated_file():
    '''
    clean  hspecies name in abdb_outfiles/abdb_segment_seq_imgt_vgene.csv
    :return:
    '''
    infile = 'abdb_outfiles/abdb_segment_seq_imgt_vgene.csv'
    df = pd.read_csv(infile)
    specs = []
    unknowns = ['nan', '-']
    for i, row in df.iterrows():
        if str(row.hspecies) not in unknowns:
            spec = '-'.join([item.strip() for item in row.hspecies.split(',')])
            if spec == 'RATTUS':
                spec = 'RATTUS RATTUS'
            elif spec == 'HOMO':
                spec = 'HOMO SAPIENS'
            elif spec == 'LAMA GLAMA':
                spec = 'LLAMA GLAMA'
        else:
            spec = 'MISSING'
        specs.append(spec)
    df['species'] = specs
    df.to_csv(infile, index=False)


def normalize_by_surface_residues():
    '''
    normalize rescount by surface residue
    :return:
    '''
    infile = 'abdb_outfiles/respairs_absort_abresnumi_segments.csv'
    sfilel = '../datasets/surface_residues/surface_residue_count_L.csv'
    sfileh = '../datasets/surface_residues/surface_residue_count_H.csv'
    df = pd.read_csv(infile)
    sdfl  = pd.read_csv(sfilel)
    sdfh  = pd.read_csv(sfileh)
    sdfldict = dict([(sdfl.iloc[i,1],sdfl.iloc[i,0]) for i in range(sdfl.shape[0])])
    sdfhdict = dict([(sdfh.iloc[i,1],sdfh.iloc[i,0]) for i in range(sdfh.shape[0])])
    print(sdfldict)
    print(sdfhdict)
    # sys.exit()
    norm_abs = []
    norm_ags  = []
    for i, row in df.iterrows():
        segment = row.segment
        sdf = df[df.segment==segment]
        antibody_counts = sdf.abres.value_counts()
        antigen_counts = sdf.agres.value_counts()
        abres = row.abres
        agres = row.agres
        chain = row.abchain
        if chain == 'L':
            norm_ab = antibody_counts[abres]/sdfldict[abres]
            norm_ag = antigen_counts[agres]/sdfldict[agres]
        elif chain == 'H':
            norm_ab = antibody_counts[abres]/sdfhdict[abres]
            norm_ag = antigen_counts[agres]/sdfhdict[agres]
        norm_abs.append(norm_ab)
        norm_ags.append(norm_ag)
    df['norm_ab'] = norm_abs
    df['norm_ag'] = norm_ags
    print(df.head())
    df.to_csv(infile, index=False)


def nlfree(l,longest_motif):
    '''
    calculate number of motif per victors formula
    :return:
    '''
    nlsum = 1
    i_s = math.floor(l/2) + l%2-1
    n= longest_motif+1
    for i in range(i_s):
        i = i+1
        if l > 2:
            a = scipy.special.comb(l-i-1,i)
        else:
            a = 1
        nl = a*(n**i)
        nlsum += nl
    return nlsum


def nlrestrict(l,longest_motif):
    '''
    calculate number of motif per victors formula
    :return:
    '''
    nlsum = 1
    i_s = math.floor(l/2) + l%2-1
    n= longest_motif+1
    n = n-l+1
    for i in range(i_s):
        i = i+1
        if l > 2:
            a = scipy.special.comb(l-i-1,i)
        else:
            a =1
        nl = a*((math.ceil(n/i))**i)
        nlsum += nl
    return nlsum

def motif_diversity(infile):
    '''
    estimate motif diversity
    :return:
    '''
    # infile = 'abdb_outfiles/respairs_paratope_segment_notationx.csv'
    df = pd.read_csv(infile)
    nlfrees = []
    nlrestricts = []
    observeds = []
    lens = []
    print(df.head())
    longest_motif = 17
    if 'epitope' in infile:
        longest_motif = 21
    for i in range(1,longest_motif):
        nlsum = nlfree(i, longest_motif)
        nlrsum = nlrestrict(i, longest_motif)
        nlfrees.append(nlsum)
        nlrestricts.append(nlrsum)
        ldf = df[df.motif_len == i]
        observed_motifs = ldf.gap_patterns.unique()
        observeds.append(len(observed_motifs))
        lens.append(i)
    diversities = nlfrees + nlrestricts + observeds
    tlens = lens*3
    limits = ['upper']*len(nlfrees) + ['lower']*len(nlrestricts) + ['observed']*len(observeds)
    print(diversities)
    print(tlens)
    print(limits)
    mdf  = pd.DataFrame({'diversity':diversities, 'limit':limits,'length':tlens})
    # print(mdf)
    # uppermdf =  mdf[mdf.limit == 'upper']
    # print(sum(uppermdf.diversity))
    # sys.exit()
    outname = infile[:-4] + '_motif_diversity.csv'
    # outname = 'abdb_outfiles/motif_diversity.csv'
    mdf.to_csv(outname, index=False)

def get_upper_overlap():
    '''
    calculate ab ag motif  upper limit overlap
    :return:
    '''
    epifile = 'abdb_outfiles/respairs_epitope_segment_notationx_motif_diversity.csv'
    parafile = 'abdb_outfiles/respairs_paratope_segment_notationx_motif_diversity.csv'
    abdf = pd.read_csv(parafile)
    abdf = abdf[abdf.limit == 'upper']
    agdf = pd.read_csv(epifile)
    agdf = agdf[agdf.limit == 'upper']
    ab_denominator = sum(abdf.diversity)
    ag_denominator = sum(agdf.diversity)
    print(ab_denominator, ag_denominator)
    join_probab = (1/ab_denominator)*(1/ag_denominator)
    n = 10**6
    njoin_probab = n*join_probab
    print(njoin_probab, join_probab)

def make_gap_dataset(infile, outdir):
    '''

    :return:
    '''
    # infile = 'abdb_outfiles/respairs_paratope_segment.csv'
    # infile = 'abdb_outfiles/respairs_epitope_segment.csv'
    df = pd.read_csv(infile)
    if 'epitope' in infile:
        df.rename({'egapset':'gapset'})
    gapdata = []
    for i, row in df.iterrows():
        try:
            gaps = row.gapset.split('-')
            for gap in gaps:
                datum = [row.pdbid, row.segment, gap, row.abchain]
                gapdata.append(datum)
        except:
            print('no gap found in %s' % row.pdbid)
    gapdf = pd.DataFrame(gapdata, columns=['pdbid', 'segment', 'gap', 'abchain'])
    outfile = infile[:-4] + '_gaponly.csv'
    inname = infile.split('/')[-1].split('.')[0]
    outfile = outdir + '/' + inname + '_gaponly.csv'
    gapdf.to_csv(outfile, index=False)

def get_levenshtein_segments():
    '''
    get levenshtein distance per segment
    :return:
    '''
    infile = 'abdb_outfiles/abdb_segment_seq_imgt_vgene.csv'
    df = pd.read_csv(infile)
    data = []
    for segment in df.segment.unique():
        segdf = df[df.segment == segment]
        print(segment)
        print(segdf.shape)
        counter = 0
        for i, row in segdf.iterrows():
            counter += 1
            print(counter)
            seq1 = row.segment_seq
            pdbid = row.pdbid
            for i2, row2 in segdf.iterrows():
                pdbid2 = row2.pdbid
                if pdbid != pdbid2:
                    seq2 = row2.segment_seq
                    ld = jellyfish.levenshtein_distance(seq1,seq2)
                    datum = [pdbid, pdbid2, segment, seq1, seq2,ld]
                    data.append(datum)
    colnames = ['pdbid1', 'pdbid2', 'segment', 'seq1', 'seq2', 'ld']
    lddf = pd.DataFrame(data, columns=colnames)
    print(lddf.head())
    outname = infile[:-4] + '_ld.csv'
    print(outname)
    lddf.to_csv(outname, index=False)

def motif_diversity_bypermutation(l):
    '''
    get possible position by permutation
    :return:
    '''
    # l = 10
    # print(list(range(3)))
    terms = math.floor(l/2) + l%2-1
    positions = list(range(1,l-1))
    coef_term = []
    for term in range(1,terms+1):
        c = list(itertools.combinations(positions,term))
        if term !=1:
            gapped_c = []
            gapped_post = []
            for posts in c:
                gaps = []
                # print(posts)
                for i,post in enumerate(posts[:-1]):
                    gap = posts[i+1] - post
                    gaps.append(gap)
                if 1 not in gaps:
                    gapped_c.append(gaps)
                    gapped_post.append(posts)
            # print(len(gapped_c), gapped_c, 'term %s, l %s'% (term,l))
            # print(gapped_post)
            coef = len(gapped_c)
            outpost = gapped_post
        else:
            # print(len(c), c, 'term %s, l %s'% (term,l))
            coef = len(c)
            outpost = c
        coef_term.append([coef, term, outpost, l])
    df = pd.DataFrame(coef_term, columns=['coef', 'term', 'positions', 'l'])
    return df

def test_motif_diversity():
    '''
    test out the algorithm
    :return:
    '''
    dfs = []
    for i in range(3,14):
        df = motif_diversity_bypermutation(i)
        dfs.append(df)
    dfs = pd.concat(dfs)
    fdata  = []
    for l in dfs.l.unique():
        ldf = dfs[dfs.l==l]
        coefs =ldf.coef
        terms = ldf.term
        print(coefs)
        formula = ''
        for coef,term in zip(coefs,terms):
            formula += '%sn^%s' % (coef,term) + '+'
        formula = formula[:-1]
        fdata.append([formula, l])
    fddf = pd.DataFrame(fdata, columns=['form', 'l'])
    print(fddf)

def naive_vgenes():
    '''
    get naive v gene distribution
    :return:
    '''
    mouse_file = 'abdb_outfiles/vgene_data_flat_list.csv'
    mdf = pd.read_csv(mouse_file).rename(columns={'Unnamed: 0': 'name'})
    mdfprebc = mdf.name.str.contains('prebc')
    mdf = mdf[mdfprebc]
    infile = 'abdb_outfiles/dewitt_plos_D1_Nb_len15.csv'
    df = pd.read_csv(infile)
    vgene_counts  = df.vGeneName.value_counts()
    # print(vgene_counts)
    abdb_infile = 'abdb_outfiles/abdb_segment_seq_imgt_vgene.csv'
    abdbdf = pd.read_csv(abdb_infile)
    habdb_vgenecounts = abdbdf[abdbdf.imgt_species=='Homo sapiens'].imgt_vgenename.value_counts()
    mabdb_vgenecounts = abdbdf[abdbdf.imgt_species=='Mus musculus_C57BL/6'].imgt_vgenename.value_counts()
    # print(abdb_vgenecounts.index)
    # print(vgene_counts.index)
    data = []
    print([item for item in vgene_counts.index if 'L' in item])
    for vgene in habdb_vgenecounts.index:
        try:
            abdb_count = habdb_vgenecounts[vgene]
            dewitt_vgene = vgene[:4] + '0' + vgene[4:]
            naive_count = vgene_counts[dewitt_vgene]
            print(abdb_count, naive_count, vgene)
            datum = [vgene, abdb_count, naive_count, 'Homo sapiens']
            data.append(datum)
        except:
            # print('%s not in naive' % vgene)
            pass
    # print(mabdb_vgenecounts.index)
    # print()
    for vgene in mabdb_vgenecounts.index:
        try:
            abdb_count = mabdb_vgenecounts[vgene]
            # naive_count = mdf[(mdf.k == vgene ) & (mdf.name.str.contains('healthy'))]#.Freq.mean()
            naive_count = mdf[(mdf.k == vgene ) & (mdf.name.str.contains('healthy'))].Freq.mean()
            # print(naive_count)
            # sys.exit()
            print(abdb_count, naive_count, vgene)
            datum = [vgene, abdb_count, naive_count, 'Mus musculus']
            data.append(datum)
        except:
            # print('%s not in naive' % vgene)
            pass
    dfout = pd.DataFrame(data, columns=['vgene', 'abdb_count', 'naive_count', 'species']).dropna()
    outname = 'abdb_outfiles/vgenes_humans_mice.csv'
    dfout.to_csv(outname, index=False)
    print(dfout[dfout.species=='Mus musculus'].corr())
    print(dfout.corr())

def shared_tenpc(infile):
    '''
    get motifs that are shared 10pc or more
    :return:
    '''
    # infile = 'abdb_outfiles/respairs_paratope_segment_notationx.csv'
    infile2 = 'abdb_outfiles_2019/respairs_absort_cutoff4_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx_count.csv'
    df = pd.read_csv(infile)
    cutoff = math.floor(df[df.segment =='CDR-H3'].raw_count.sum()/10)
    shareddf = df[df.raw_count >= cutoff]
    outname = infile[:-4] + '_shared.csv'
    shareddf.to_csv(outname, index=False)


def get_angle():
    '''
    gets angels from paratope/epitope
    :return:
    '''
    infile = 'abdb_outfiles/respairs_segment_notationx_merged.csv'
    df = pd.read_csv(infile)
    # df = df[df.pdbid == '1IKF_1']
    # print(df.head())
    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    aadictr = dict([(single,triple) for triple, single in aadict.items()])
    data = []
    for pdbid in df.pdbid.unique():
        pdbdf = df[df.pdbid == pdbid]
        print(pdbdf)
        pdbfile = '../datasets/NR_LH_Protein_Martin/' + pdbid + '.pdb'
        pdbcontents = open(pdbfile).read().splitlines()
        pdbdict = {}
        pdbresdict = {}
        for content in pdbcontents:
            if content.startswith('ATOM'):
                resname = content[17:20].strip()
                chainid = content[21].strip()
                resseq = content[22:26].strip() + content[26].strip()
                x,y,z = float(content[30:38]), float(content[38:46]), float(content[46:54])
                coords = [x,y,z]
                key = resname + chainid + resseq
                if key not in pdbdict:
                    pdbdict[key] = [coords]
                    pdbresdict[key] = content + '\n'
                else:
                    pdbdict[key].append(coords)
                    pdbresdict[key] += content + '\n'
        for i, row in pdbdf.iterrows():
            # print(row)
            plen = row.plen
            elen = row.elen
            print(plen)
            if plen >= 3:
                print('hey')
                paratope = list(row.paratope)
                abresnumis = row.abresnumiset.split('-')
                chains = [row.abchain]*len(paratope)
                centers = []
                pdbres= ''
                pseudos = ''
                for res,chain,abresnumi in zip(paratope, chains, abresnumis):
                    key = ''.join([aadictr[res],chain,abresnumi])
                    coordinates = np.array(pdbdict[key])
                    center = np.sum(coordinates, axis=0)/len(coordinates)
                    centers.append(center)
                    pdbres += pdbresdict[key]
                start, mid, end,  = centers[0], centers[1:-1], centers[-1]
                mid = np.sum(mid, axis=0)/len(mid)
                startcoord = ['{:8.4f}'.format(item) for item in start]
                midcoord = ['{:8.4f}'.format(item) for item in mid]
                endcoord = ['{:8.4f}'.format(item) for item in end]
                pseudo_start = 'ATOM'.ljust(6) + content[6:17] + 'PSE' + content[20:30] + startcoord[0] + startcoord[1]+ \
                               startcoord[2] + content[54:]
                pseudo_mid = 'ATOM'.ljust(6) + content[6:17] + 'PSE' + content[20:30] + midcoord[0] + midcoord[1]+ \
                               midcoord[2] + content[54:]
                pseudo_end = 'ATOM'.ljust(6) + content[6:17] + 'PSE' + content[20:30] + endcoord[0] + endcoord[1]+ \
                             endcoord[2] + content[54:]
                pseudos = '\n'.join([pseudo_start,pseudo_mid, pseudo_end])
                v1 = mid-start
                v2 = mid-end
                costheta = np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
                arccostheta = math.degrees(math.acos(costheta))
                pdbresoutname = 'abdb_outfiles/interacting_pdb/' + '%s_%s.pdb' % (pdbid, row.segment)
                pseudooutname = 'abdb_outfiles/interacting_pdb/' + '%s_%s_pseudo.pdb' % (pdbid, row.segment)
                pdbresoutfile = open(pdbresoutname, 'w')
                pseudooutfile = open(pseudooutname, 'w')
                pdbresoutfile.write(pdbres)
                pseudooutfile.write(pseudos)
                p_angle = round(arccostheta)
            elif plen <3:
                p_angle = np.nan
            if elen >= 3:
                paratope = list(row.epitope)
                abresnumis = row.agresnumiset.split('-')
                chains = [row.agchain]*len(paratope)
                centers = []
                pdbres= ''
                pseudos = ''
                for res,chain,abresnumi in zip(paratope, chains, abresnumis):
                    key = ''.join([aadictr[res],chain,abresnumi])
                    coordinates = np.array(pdbdict[key])
                    center = np.sum(coordinates, axis=0)/len(coordinates)
                    centers.append(center)
                    pdbres += pdbresdict[key]
                start, mid, end,  = centers[0], centers[1:-1], centers[-1]
                mid = np.sum(mid, axis=0)/len(mid)
                startcoord = ['{:8.4f}'.format(item) for item in start]
                midcoord = ['{:8.4f}'.format(item) for item in mid]
                endcoord = ['{:8.4f}'.format(item) for item in end]
                pseudo_start = 'ATOM'.ljust(6) + content[6:17] + 'PSE' + content[20:30] + startcoord[0] + startcoord[1]+ \
                               startcoord[2] + content[54:]
                pseudo_mid = 'ATOM'.ljust(6) + content[6:17] + 'PSE' + content[20:30] + midcoord[0] + midcoord[1]+ \
                             midcoord[2] + content[54:]
                pseudo_end = 'ATOM'.ljust(6) + content[6:17] + 'PSE' + content[20:30] + endcoord[0] + endcoord[1]+ \
                             endcoord[2] + content[54:]
                pseudos = '\n'.join([pseudo_start,pseudo_mid, pseudo_end])
                v1 = mid-start
                v2 = mid-end
                costheta = np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
                arccostheta = math.degrees(math.acos(costheta))
                pdbresoutname = 'abdb_outfiles/interacting_pdb/' + '%s_%sag.pdb' % (pdbid, row.segment)
                pseudooutname = 'abdb_outfiles/interacting_pdb/' + '%s_%s_pseudoag.pdb' % (pdbid, row.segment)
                pdbresoutfile = open(pdbresoutname, 'w')
                pseudooutfile = open(pseudooutname, 'w')
                pdbresoutfile.write(pdbres)
                pseudooutfile.write(pseudos)
                e_angle = round(arccostheta)
            elif elen <3:
                e_angle = np.nan
            datum = row.tolist() + [p_angle, e_angle]
            data.append(datum)
            print(datum)
    columns = row.index.tolist() + ['p_angle', 'e_angle']
    outdf = pd.DataFrame(data, columns=columns)
    outname = infile[:-4] + '_angle.csv'
    outdf.to_csv(outname,index=False)

def check_angle():
    '''

    :return:
    '''
    infile = 'abdb_outfiles/respairs_segment_notationx_merged_angle.csv'
    df = pd.read_csv(infile)
    print(df.e_angle.dropna().shape)
    print(df[df.elen >=3].shape)

def prepdata_network():
    '''
    prep data for network analysis
    :return:
    '''
    infile = 'abdb_outfiles/respairs_segment_notationx_merged_angle.csv'
    df = pd.read_csv(infile)
    print(df.info())
    motif_pairs = {}
    for i, row in df.iterrows():
        abmotif = row.ab_gap_pattern
        agmotif = row.ag_gap_pattern
        motif_pair = '-'.join([abmotif,agmotif])
        if motif_pair not in motif_pairs:
            motif_pairs[motif_pair] = 1
        else:
            motif_pairs[motif_pair] += 1
    data = []
    for pair, value in motif_pairs.items():
        ab, ag  = pair.split('-')
        datum = [ab, ag, value]
        data.append(datum)
    df = pd.DataFrame(data, columns = ['abmotif', 'agmotif', 'count'])
    outname = 'abdb_outfiles/motifpairs.csv'
    df.to_csv(outname, index=False)


def prep_meta_heatmap():
    '''
    prep meta data for heatmap
    :return:
    '''
    infile2 = 'abdb_outfiles/respairs_segment_notationx_merged_angle.csv'
    df = pd.read_csv(infile2)
    dfab = df.drop_duplicates(subset='ab_gap_pattern')
    ab_gap_statuses = []
    ag_gap_statuses = []
    print(dfab.head())
    for i, row in dfab.iterrows():
        if row.ab_gapstatus ==0:
            ab_gap_statuses.append('Continuous')
        else:
            ab_gap_statuses.append('Discontinuous')
    dfag = df.drop_duplicates(subset='ag_gap_pattern')
    for i, row in dfag.iterrows():
        if row.ag_gapstatus ==0:
            ag_gap_statuses.append('Continuous')
        else:
            ag_gap_statuses.append('Discontinuous')
    abmeta = pd.DataFrame({'Gap Status':ab_gap_statuses}).set_index(dfab.ab_gap_pattern)
    agmeta = pd.DataFrame({'Gap Status':ag_gap_statuses}).set_index(dfag.ag_gap_pattern)
    abmeta.to_csv('abdb_outfiles/abmeta.csv')
    agmeta.to_csv('abdb_outfiles/agmeta.csv')

def prepdata_abag_residue_distribution():
    '''

    :return:
    '''
    infile = 'abdb_outfiles/respairs_absort_abresnumi_segments.csv'
    df = pd.read_csv(infile)
    print(df.head())
    abdfs = []
    agdfs = []
    for segment in df.segment.unique():
        print(segment)
        sdf = df[df.segment == segment]
        sdf = sdf.drop_duplicates(subset=['abres'])
        agsdf = sdf.drop_duplicates(subset='agres')
        abdfs.append(sdf)
        agdfs.append(agsdf)
    abdf = pd.concat(abdfs)
    agdf = pd.concat(agdfs)
    outname1 =  infile.split('.')[0] + '_residue_freq_ab.csv'
    outname2 =  infile.split('.')[0] + '_residue_freq_ag.csv'
    abdf.to_csv(outname1)
    agdf.to_csv(outname2)
    sys.exit()
    print(abdf.head())

def random_motif_overlap():
    '''
    get random motif overlap
    :return:
    '''
    infile_ab = 'abdb_outfiles/respairs_paratope_segment_notationx_top-1_gap_patterns.csv'
    infile_ag = 'abdb_outfiles/respairs_epitope_segment_notationx_top-1_gap_patterns.csv'
    dfab = pd.read_csv(infile_ab)
    dfag = pd.read_csv(infile_ag)
    abmotif = dfab.gap_patterns.tolist()
    agmotif = dfag.gap_patterns.tolist()
    overlaps = []
    for i in range(100):
        rab = random.sample(set(abmotif), 100)
        rag = random.sample(set(agmotif), 370)
        overlap = list(set(rab) & set(rag))
        overlap_fraction  = len(overlap)/len(set(rab))
        # print(len(set(rab))/len(set(rag)))
        overlaps.append(overlap_fraction)
    print(sum(overlaps)/ len(overlaps))

def check_dynet_data():
    '''
    checkout dynet data
    :return:
    '''
    infile = 'abdb_outfiles/DyNet_Central_Reference_Network_default_node.csv'
    df = pd.read_csv(infile)
    print(df.head())

def prepdata_rewiring():
    '''
    prepare rewiring dataframe
    :return:
    '''
    infiles = fifi.find_files('abdb_outfiles', 'rewiring.csv')
    infiles = [item for item in infiles if 'paired' not in item]
    data = []
    for infile in infiles:
        df = pd.read_csv(infile)
        for i, row in df.iterrows():
            residue = row['name']
            rewscore = round(row['DyNet Rewiring (Dn-score)'],2)
            # rewscore = round(row['Dn-Score (degree corrected)'],2)
            pair = '_'.join(infile.split('/')[1].split('_')[0:3])
            source = pair.split('_')[0]
            datum  = [residue, rewscore, pair, source]
            data.append(datum)
    colnames = ['residue', 'rewiring_score', 'pair', 'source']
    outdf = pd.DataFrame(data, columns=colnames)
    outname = 'abdb_outfiles/paired_motif_rewiring.csv'
    outdf.to_csv(outname, index=False)

def get_residue_from_motif(sequence, motif):
    '''
    get residue decomposition for a given motif
    :return:
    '''
    gap_indices = [i for i in range(len(motif)) if motif[i] != 'X']
    gap_sizes = [int(item) for item in motif if item != 'X']
    gap_sizes_adj = [0] + [int(item)-1 for item in motif if item != 'X']
    res_indices = [i for i in range(len(motif)) if motif[i] == 'X']
    dseqs = []
    indices = []
    if len(gap_sizes) > 0:
        newris = []
        for ri in res_indices[:]:
            new_i = ri + sum(gap_sizes_adj[:ri])
            newris.append(new_i)
        res_indices = newris
    window_len = len(res_indices) + sum(gap_sizes)
    window_indices = range(len(sequence)-window_len + 1)
    for wi in window_indices:
        segment = sequence[wi:wi+window_len]
        residues = ''.join([segment[i] for i in res_indices])
        print(segment, sequence, residues)
        print(wi)
        dseqs.append(residues)
        indices.append(wi)
    return dseqs, indices

def prepdata_dewitt():
    '''
    prep dewitt data for seq net and rewiring
    :return:
    '''
    infiles = fifi.find_files('../datasets/dewitt_2016', '.csv')
    motif = 'XXX'
    # infiles = [item for item in infiles if 'D1' not in item]
    for infile in infiles[:]:
        print(infile)
        df = pd.read_csv(infile)
        sequences = df.aminoAcid
        dseqss = []
        indicess = []
        for sequence in sequences:
            dseqs, indices = get_residue_from_motif(sequence, motif)
            dseqss += dseqs
            indicess += indices
        print(len(dseqss))
        gap_patterns = [motif]*len(dseqss)
        outdf = pd.DataFrame(dict([('paratope', dseqss), ('resindex', indicess), ('gap_patterns',gap_patterns)]))
        print(outdf.head())
        outname = infile.split('/')[-1].split('.')[0]
        print(outname)
        outfile = 'abdb_outfiles/paratope_decomposition_%s_%s.csv' % (outname, motif)
        outdf.to_csv(outfile, index=False)

def loop_get_context_pattern_top3():
    '''
    get seqnet data for dewitt stuff
    now with position
    :return:
    '''
    infiles = fifi.find_files('abdb_outfiles', 'paratope_decomposition')
    infiles = [item for item in infiles if 'next' not in item][-3:]
    # print(infiles)
    # sys.exit()
    for infile in infiles[:]:
        print(infile)
        get_context_pattern_top3_position(infile)
        # sys.exit()

def rename_dewitt_files():
    '''
    assign shorter names to the files
    :return:
    '''
    infiles = fifi.find_files('abdb_outfiles', 'len15_XXX_XXX_top3motif_respairs')
    for infile in infiles:
        if 'position.csv' in infile:
            print(infile)
            parts = infile.split('/')
            nameparts = parts[1].split('_')
            newname = parts[0] + '/' + '_'.join([nameparts[4], nameparts[5], nameparts[7], nameparts[-2],
                                                 nameparts[-1]])
            print(newname)
        parts = infile.split('/')
        nameparts = parts[1].split('_')
        newname = parts[0] + '/' + '_'.join([nameparts[4], nameparts[5], nameparts[7], nameparts[-1]])
        os.system('cp %s %s' % (infile, newname))
        print(infile)


def prepdata_rewiring_dewitt():
    '''
    prep rewiring data from dewitt seqnet
    :return:
    '''
    infiles = fifi.find_files('abdb_outfiles', 'rewiring.csv')
    infiles = [item for item in infiles if 'memory' in item or 'naive' in item]
    print(infiles)
    data = []
    for infile in infiles[:]:
        df = pd.read_csv(infile)
        for i, row in df.iterrows():
            residue = row['name']
            rewscore = round(row['DyNet Rewiring (Dn-score)'],2)
            # rewscore = round(row['Dn-Score (degree corrected)'],2)
            source = '_'.join(infile.split('/')[1].split('_')[0:-1])
            datum  = [residue, rewscore, source]
            print(datum)
            data.append(datum)
    colnames = ['residue', 'rewscore', 'source']
    outdf = pd.DataFrame(data, columns=colnames)
    outname = 'abdb_outfiles/dewitt_allrep_rewiring.csv'
    print(outdf.head())
    outdf.to_csv(outname, index=False)


def normalize_edge_pairs():
    '''
    normalize edges, rebuilt network see if this changes stuff
    :return:
    '''
    infiles = fifi.find_files('abdb_outfiles', 'XXX_next.csv') + fifi.find_files('abdb_outfiles',
                                                                                 'XXX_position.csv')
    for infile in infiles:
        df = pd.read_csv(infile)
        counts = df['count']
        print(counts)
        normcount = (counts- min(counts))/(max(counts)-min(counts))
        df['normcount'] = normcount
        # print(normcount)
        print(df.head())
        outname = infile.split('.')[0] + '_norm.csv'
        df.to_csv(outname, index=False)
        # sys.exit()


def normalize_edge_pairs_composite():
    '''
    normalize edges, rebuilt network see if this changes stuff
    for composite nets: XXX, X1X, XX
    :return:
    '''
    # infiles = fifi.find_files('abdb_outfiles', 'XXX_next.csv')
    infiles = ['abdb_outfiles/paratope_XXX_XX_X1X_top3motif_respairs_next.csv', 'abdb_outfiles/epitope_XX_X1X_XXX_top3motif_respairs_next.csv']
    for infile in infiles:
        df = pd.read_csv(infile)
        counts = df['count']
        print(counts)
        normcount = (counts- min(counts))/(max(counts)-min(counts))
        df['normcount'] = normcount
        # print(normcount)
        print(df.head())
        outname = infile.split('.')[0] + '_norm.csv'
        df.to_csv(outname, index=False)
        # sys.exit()

def merge_repertoire_edges():
    '''
    merge edges from repertoire, use position files
    :return:
    '''
    infiles = fifi.find_files('abdb_outfiles', 'position_norm.csv')
    dfs = []
    for infile in infiles:
        df = pd.read_csv(infile)
        print(df.head())
        parts = infile.split('/')[1].split('_')
        donor = parts[0]
        rep =  parts[1]
        motif = parts[2]
        df['donor'] = [donor]*df.shape[0]
        df['reptype'] = [rep]*df.shape[0]
        df['motif'] = [motif]*df.shape[0]
        df['idx'] = [item[1:] for item in df.source]
        dfs.append(df[-20:])
    outdf = pd.concat(dfs)
    outfile = 'abdb_outfiles/merged_repertoire_edges.csv'
    outdf.to_csv(outfile, index=False)

def prepdata_parepi_naivememory_rewiring():
    '''
    prepdata for paratope epitipe vs memory naive rewiring
    :return:
    '''
    infiles = fifi.find_files('abdb_outfiles', 'rewiring.csv')
    infiles = [item for item in infiles if '_memory' in item or '_naive' in item][:-1]
    data = []
    for infile in infiles:
        df = pd.read_csv(infile)
        for i, row in df.iterrows():
            residue = row['name']
            rewscore = round(row['DyNet Rewiring (Dn-score)'],2)
            # rewscore = round(row['Dn-Score (degree corrected)'],2)
            parts = infile.split('/')[-1].split('_')
            pair = '_'.join(parts[:-1])
            datum  = [residue, rewscore, pair]
            data.append(datum)
    colnames = ['residue', 'rewscore', 'pair']
    outdf  = pd.DataFrame(data, columns=colnames)
    outname =  'abdb_outfiles/paraepi_naivememory.csv'
    outdf.to_csv(outname, index=False)
    print(outdf.head())


def prepdata_decomposition_probability():
    '''
    prep data for decomposition probability plot. for XXX
    :return:
    '''
    infiles = fifi.find_files('abdb_outfiles', 'decomposition')
    infiles = [item for item in infiles if 'plos' in item and 'next' not in item]
    infiles = infiles + ['abdb_outfiles/respairs_paratope_segment_notationx.csv',
                         'abdb_outfiles/respairs_epitope_segment_notationx.csv']
    # print(infiles)
    # sys.exit()
    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    aalist= aadict.values()
    dfs= []
    dfsfull = []
    for infile in infiles:
        df = pd.read_csv(infile)[:]
        df = df[df.gap_patterns == 'XXX']
        donor = infile.split('/')[-1].split('_')[-4]
        reptype = infile.split('/')[-1].split('_')[-3]
        if 'respairs_paratope' in infile:
            donor = 'paratope'
            reptype = 'abdb'
            df = df[df.gap_patterns == 'XXX']
            print(df.head())
        elif 'respairs_epitope' in infile:
            donor = 'epitope'
            reptype = 'abdb'
            df = df[df.gap_patterns == 'XXX']
            df = df.rename(columns={'epitope':'paratope'})
        p1 = dict([(item,0) for item in aalist])
        p2 = dict([(item,0) for item in aalist])
        p3 = dict([(item,0) for item in aalist])
        posdicts = [p1,p2,p3]
        for i, row in df.iterrows():
            seqs = list(row.paratope)
            for char, posdict in zip(seqs, posdicts):
                posdict[char] += 1
        for posdict in posdicts:
            totalres = df.shape[0]
            for char in posdict:
                posdict[char] = posdict[char]/float(totalres)
        seqprobas = []
        dfunique = df.drop_duplicates(subset=['paratope'])
        print(dfunique.shape)
        print(df.shape)
        for i,row in dfunique.iterrows():
            chars = list(row.paratope)
            joinproba = 1
            for i2, char in enumerate(chars):
                joinproba*=posdicts[i2][char]
            seqprobas.append(joinproba)
        dfunique['joinprobab'] = seqprobas
        dfunique['source'] = [donor]*dfunique.shape[0]
        dfunique['reptype'] = [reptype]*dfunique.shape[0]
        dfunique = dfunique.sort_values(by='joinprobab')[-5:]
        dfuniquefull = dfunique.sort_values(by='joinprobab')
        dfs.append(dfunique)
        dfsfull.append(dfuniquefull)
    outdf = pd.concat(dfs)
    outdffull = pd.concat(dfsfull)
    outname = 'abdb_outfiles/decomposition_probability.csv'
    outnamefull = 'abdb_outfiles/decomposition_probability_full.csv'
    outdf.to_csv(outname, index=False)
    outdffull.to_csv(outnamefull, index=False)

def annotate_interaction_network_degree():
    '''
    annotate abdb_outfiles/ab_ag_network_node.csv
    :return:
    '''
    infile = 'abdb_outfiles/ab_ag_network_node.csv'
    df = pd.read_csv(infile)
    sources = []
    for i, row in df.iterrows():
        if row.SUID < 1326:
            sources.append('paratope')
        else:
            sources.append('epitope')
    df['source']  = sources
    print(df.head())
    outname =  infile.split('.')[0] + '_paraepi.csv'
    df.to_csv(outname, index=False)
    # print(df.head())

def filter_dewitt_decomposition():
    '''
    filter the first three chars out ( CAR and co)
    :return:
    '''
    infiles = fifi.find_files('abdb_outfiles', 'decomposition_dewitt_plos')
    infiles = [item for item in infiles if 'next' not in item]
    for infile in infiles:
        df = pd.read_csv(infile)
        df = df[df.resindex >2]
        print(df.head())
        df.to_csv(infile, index=False)
        print(infiles)

def make_aadict():
    '''
    makes aadict for three letters to one letter conversion
    :return:
    '''

    aafile = '../datasets/amino_acids/the_twenty.txt'
    aacontent = open(aafile).read().splitlines()
    aadict = dict((item.split()[1].upper(),item.split()[-1]) for item in aacontent)
    return aadict

def prepdata_dlbit():
    '''
    preps the data for deeplearning bit. looks only on CDR-H3
    :return:
    '''
    motif = 'XXX'
    segment = 'CDR-H3'
    infile = 'abdb_outfiles/respairs_paratope_segment_notationx.csv'
    df = pd.read_csv(infile)
    # df = df[df.gap_patterns == motif]
    df = df[df.segment == segment]
    # print(df.gap_patterns.value_counts())
    print(df.gap_patterns.value_counts())
    print(df.shape)
    strseqs  = []
    strresseqs = []
    lenstrseqs = []
    aadict = make_aadict()
    for i,row in df.iterrows():
        pdbid = row.pdbid
        print(pdbid)
        pdbfile = '../datasets/NR_LH_Protein_Martin/' + pdbid + '.pdb'
        contents = open(pdbfile).readlines()
        # get CDR-H3 residues, chain H, resnum range 95-102
        seqs = []
        resseqs = []
        for line in contents:
            if line.startswith('ATOM'):
                chain = line[21]
                resseq = int(line[22:26])
                if chain == 'H' and 95<=resseq<=102:
                    res = line[17:20]
                    aa = aadict[res]
                    resseq_achar = str(resseq) + line[26].strip()
                    if resseq_achar not in resseqs:
                        seqs.append(aa)
                        resseqs.append(resseq_achar)
        strseq = ''.join(seqs)
        strresseq = '-'.join(resseqs)
        lenstrseq = len(strseq)
        strseqs.append(strseq)
        strresseqs.append(strresseq)
        lenstrseqs.append(lenstrseq)

    df['seq'] = strseqs
    df['resseq'] =  strresseqs
    df['seqlen'] = lenstrseqs
    outfile  = '../dl/dataset/abdb_cdrh3.csv'
    df.to_csv(outfile, index=False)
    print(df.head())

def get_motif_len():
    '''
    gets motif lens for the merged file
    :return:
    '''
    infile = 'abdb_outfiles/respairs_segment_notationx_merged_angle.csv'
    df = pd.read_csv(infile).iloc[:]
    abmotiflens = []
    agmotiflens = []
    abmotiflengaps = []
    agmotiflengaps = []
    for i, row in df.iterrows():
        abmotif = row.ab_gap_pattern
        abparts = abmotif.split('X')
        agmotif = row.ag_gap_pattern
        agparts = agmotif.split('X')
        # print(agmotif, agparts)
        if row.ab_gapstatus != 0:
            abmotiflen = len(abparts)
            abints = [int(item) for item in abparts if item != '']
            abmotiflengap = abmotiflen -1 + sum(abints)
        else:
            abmotiflen = len(abmotif)
            abmotiflengap = len(abmotif)
        if row.ag_gapstatus !=0:
            agmotiflen = len(agparts)
            agints = [int(item) for item in agparts if item != '']
            agmotiflengap = agmotiflen -1 + sum(agints)
        else:
            agmotiflen = len(agmotif)
            agmotiflengap = len(agmotif)
        abmotiflens.append(abmotiflen)
        agmotiflens.append(agmotiflen)
        abmotiflengaps.append(abmotiflengap)
        agmotiflengaps.append(agmotiflengap)
    df['ab_motiflen'] = abmotiflens
    df['ag_motiflen'] = agmotiflens
    df['ab_motiflengap'] = abmotiflengaps
    df['ag_motiflengap'] = agmotiflengaps
    #remove outliers
    df = df[df.ag_motiflengap < 1000]
    outfile = infile.split('.')[0] + '_len.csv'
    df.to_csv(outfile, index=False)


def manuscript_graphics():
    '''
    plots for manucscripts
    :return:
    '''
    # plot_residue_number_distribution('abdb_outfiles/respairs_absort_abresnumi_segments.csv',
    # 								 'interacting_residue_position')
    # plot_interacting_residue_distribution('abdb_outfiles/respairs_absort_abresnumi_segments.csv',
    # 									  'interacting_residue_distribution')

    # plot_interacting_residue_distribution_surface('abdb_outfiles/respairs_absort_abresnumi_segments.csv',
    # 									  'interacting_residue_distribution_surface_normalized')
    # plot_paratope_length('abdb_outfiles/respairs_paratope_segment.csv', 'paratope_length_segment_wise')
    # plot_paratope_gapsize('abdb_outfiles/respairs_paratope_segment.csv', 'paratope_gapsize_segment_wise')
    # plot_gap_patterns_gapset('abdb_outfiles/respairs_paratope_segment.csv', 'paratope_gap_patterns_gapset')
    # plot_gap_patterns_position('abdb_outfiles/respairs_paratope_segment_position_residue.csv',
    #                            'paratope_gap_patterns_position')
    # plot_gap_patterns_position_chain('abdb_outfiles/respairs_paratope_segment_position_residue.csv',
    #                            'paratope_gap_patterns_position')
    # add_gap_patterns_position_residue('abdb_outfiles/respairs_paratope_segment.csv',
    # 								  'abdb_outfiles/respairs_paratope_segment_position_residue.csv')
    # plot_residue_position_cdr3('abdb_outfiles/respairs_paratope_segment.csv')
    # plot_gap_patterns_distribution('abdb_outfiles/respairs_paratope_segment.csv')
    # plot_gap_pos_patterns_distribution('abdb_outfiles/respairs_paratope_segment.csv')
    # pdb2seqtable()
    # add_notationx('abdb_outfiles/respairs_paratope_segment.csv', 'gapset', 'abresnumiset')
    # add_notationx('abdb_outfiles/respairs_epitope_segment.csv', 'egapset', 'agresnumiset')
    # plot_gap_patterns_all_segments('abdb_outfiles/respairs_paratope_segment_notationx.csv')
    # plot_gap_pos_patterns_all_segments('abdb_outfiles/respairs_paratope_segment_notationx.csv')
    # get_pattern_residue_decomposition('abdb_outfiles/respairs_paratope_segment_notationx.csv', 'paratope',
    #                                   'abresnumiset', 'CDR-H3')
    # prepdata_ag()
    # add_notationx('abdb_outfiles/respairs_epitope_segment.csv', 'egapset', 'agresnumiset')
    # plot_gap_patterns_all_segments('abdb_outfiles/respairs_epitope_segment_notationx.csv')
    # get_pattern_residue_decomposition('abdb_outfiles/respairs_epitope_segment_notationx.csv', 'epitope',
    #                                   'agresnumiset', 'CDR-H3')
    # plot_epitope_length('abdb_outfiles/respairs_epitope_segment_notationx.csv')
    # plot_epitope_gapsize('abdb_outfiles/respairs_epitope_segment_notationx.csv')
    # plot_pattern_residue_decomposition_cdrs_ab()
    # plot_pattern_residue_decomposition_cdrs_ag()
    # prepdata_pattern_interaction_map('abdb_outfiles/respairs_paratope_segment_notationx.csv',
    #                                  'abdb_outfiles/respairs_epitope_segment_notationx.csv')
    # plot_pattern_interaction_map('abdb_outfiles/respairs_segment_notationx_merged.csv', 'all')
    # plot_pattern_interaction_map_all_segments()
    # plot_residue_interaction_map('abdb_outfiles/respairs_absort_abresnumi_segments.csv')
    # plot_residue_position_interaction_map('abdb_outfiles/respairs_absort_abresnumi_segments.csv', 'all', 'H')
    # plot_gap_pattern_pwm('abdb_outfiles/respairs_paratope_segment_notationx.csv', 'XXXXXX', 'CDR-H3')
    # loop_plot_gap_pattern_pwm()
    # get_context_overview()
    # get_context_pattern('CDR-H3')

def manuscript_graphics2():
    '''
    plots graphics shown in the gdoc
    :return:
    '''
    # plot_gap_pattern_pwm('abdb_outfiles/respairs_paratope_segment_notationx.csv', 'XXXXXX', 'CDR-H3')
    # plot_gap_pattern_pwm('abdb_outfiles/respairs_epitope_segment_notationx.csv', 'XX', 'CDR-H3')
    # get_context_pattern('abdb_outfiles/respairs_paratope_segment_notationx.csv','CDR-H3')
    # get_context_pattern('abdb_outfiles/respairs_epitope_segment_notationx.csv','CDR-H3')
    # plot_residue_interaction_map('abdb_outfiles/respairs_absort_abresnumi_segments.csv', 'CDR-H3')
    # plot_residue_position_interaction_map('abdb_outfiles/respairs_absort_abresnumi_segments.csv', 'all', 'H')
    # get_ab_ag_motif_pairs('abdb_outfiles/respairs_segment_notationx_merged.csv')
    # plot_gap_patterns_all_segments('abdb_outfiles/respairs_paratope_segment_notationx.csv')
    # plot_gap_patterns_all_segments('abdb_outfiles/respairs_epitope_segment_notationx.csv')
    # plot_pattern_interaction_map_all_segments()
    # merge_segment_motif('abdb_outfiles/respairs_paratope_segment_notationx.csv')
    # get_levenshtein_distance('abdb_outfiles/respairs_paratope_segment_notationx_motif_merged.csv')
    # add_motif_stats('abdb_outfiles/respairs_paratope_segment_notationx.csv')
    # add_motif_stats('abdb_outfiles/respairs_epitope_segment_notationx.csv')
    # get_full_segment_seq('abdb_outfiles/abdb_segment.csv')
    # get_vgenes()
    # annotate_vgenes()
    # clean_annotated_file()
    # normalize_by_surface_residues()
    # make_gap_dataset()
    # motif_diversity('abdb_outfiles/respairs_paratope_segment_notationx.csv')
    # motif_diversity('abdb_outfiles/respairs_epitope_segment_notationx.csv')
    # get_levenshtein_distance('abdb_outfiles/respairs_paratope_segment_notationx_motif_merged.csv')
    # get_levenshtein_segments()
    # motif_diversity_bypermutation()
    # test_motif_diversity()
    # naive_vgenes()
    # shared_tenpc('abdb_outfiles/respairs_paratope_segment_notationx_top-1_gap_patterns.csv')
    # shared_tenpc('abdb_outfiles/respairs_epitope_segment_notationx_top-1_gap_patterns.csv')
    # get_angle()
    # check_angle()
    # prepdata_network()
    # prep_meta_heatmap()
    # get_context_pattern_top3('abdb_outfiles/respairs_paratope_segment_notationx.csv')
    # get_context_pattern_top3('abdb_outfiles/respairs_epitope_segment_notationx.csv')
    # prepdata_abag_residue_distribution()
    # random_motif_overlap()
    # check_dynet_data()
    # prepdata_rewiring()
    # prepdata_dewitt()
    # get_context_pattern_top3('abdb_outfiles/paratope_dewitt_plos_D1_Nb_len15_XXX.csv')
    # loop_get_context_pattern_top3()
    # rename_dewitt_files()
    # prepdata_rewiring_dewitt()
    # normalize_edge_pairs()
    # merge_repertoire_edges()
    # normalize_edge_pairs_composite()
    # prepdata_parepi_naivememory_rewiring()
    # prepdata_decomposition_probability()
    # annotate_interaction_network_degree()
    # filter_dewitt_decomposition()
    # prepdata_dlbit()
    # get_motif_len()
    # get_upper_overlap()


#grid arrange
#run stuff
#prepdata_ab()
#plotdata()
#check_data()
# manuscript_graphics()
manuscript_graphics2()
