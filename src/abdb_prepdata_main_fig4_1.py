# import stuff
import os
import sys
import pandas as pd


def cross_reactivity_paratope():
    '''
    prep neat data for cross reactivity plots
    :return:
    '''
    peinfile = 'abdb_outfiles_2019/paratope_epitope_internet_edges.csv'
    rinfile = 'abdb_outfiles_2019/random_internet_edges.csv'
    peinfile_nodes = 'abdb_outfiles_2019/paratope_epitope_internet_nodes.csv'
    nodesdf = pd.read_csv(peinfile_nodes)
    pedf = pd.read_csv(peinfile).iloc[:]
    rdf = pd.read_csv(rinfile)
    print(pedf.head())
    data = []
    nodes = nodesdf.id
    for motif in nodes:
        mdf = pedf[pedf.source == motif]
        partners = mdf.target
        datum  = [motif]
        for motif2 in nodes:
            mdf2 = pedf[pedf.source == motif2]
            partners2 = mdf2.target
            intersect = set(partners) & set(partners2)
            percent_overlap = round(len(intersect)/float(len(partners))*100, 1)
            # print(motif, percent_overlap)
            datum.append(percent_overlap)
        data.append(datum)
    colnames = ['motif'] + pedf.source.unique().tolist()
    print(colnames)
    outdf = pd.DataFrame(data, columns=colnames)
    outname = peinfile.split('.')[0] + '_paratope_cross.csv'
    print(outname)
    outdf.to_csv(outname, index=False)


def cross_reactivity_paratope_topn():
    '''
    prep neat data for cross reactivity plots
    :return:
    '''
    peinfile = 'abdb_outfiles_2019/paratope_epitope_internet_edges.csv'
    peinfile_nodes = 'abdb_outfiles_2019/paratope_epitope_internet_nodes.csv'
    pedf = pd.read_csv(peinfile).iloc[:]
    nodesdf = pd.read_csv(peinfile_nodes)
    print(nodesdf.head())
    print(pedf.head())
    data = []
    n = 10
    nodes = nodesdf.id[:n]
    for motif in nodes:
        mdf = pedf[pedf.source == motif]
        partners = mdf.target
        datum  = [motif]
        for motif2 in nodes:
            mdf2 = pedf[pedf.source == motif2]
            partners2 = mdf2.target
            intersect = set(partners) & set(partners2)
            percent_overlap = round(len(intersect)/float(len(partners))*100, 1)
            # print(motif, percent_overlap)
            datum.append(percent_overlap)
        data.append(datum)
    colnames = ['motif'] + nodes.tolist()
    print(colnames)
    outdf = pd.DataFrame(data, columns=colnames)
    outname = peinfile.split('.')[0] + '_paratope_cross_top%s.csv' % n
    print(outname)
    outdf.to_csv(outname, index=False)


def cross_reactivity_density():
    '''
    prep neat data for cross reactivity density plots both paratope and epitope
    :return:
    '''
    peinfile = 'abdb_outfiles_2019/paratope_epitope_internet_edges.csv'
    rinfile = 'abdb_outfiles_2019/random_internet_edges.csv'
    peinfile_nodes = 'abdb_outfiles_2019/paratope_epitope_internet_nodes.csv'
    rinfile_nodes = 'abdb_outfiles_2019/random_internet_nodes.csv'
    penodesdf = pd.read_csv(peinfile_nodes)
    rnodesdf = pd.read_csv(rinfile_nodes)
    pedf = pd.read_csv(peinfile).iloc[:]
    rdf = pd.read_csv(rinfile)
    print(pedf.head())
    data = []
    nodes_paratope = [item for item in penodesdf.id.tolist() if '*' not in item][:]
    for motif in nodes_paratope:
        mdf = pedf[pedf.source == motif]
        partners = mdf.target
        for motif2 in nodes_paratope:
            mdf2 = pedf[pedf.source == motif2]
            partners2 = mdf2.target
            intersect = set(partners) & set(partners2)
            percent_overlap = round(len(intersect)/float(len(partners))*100, 1)
            # print(motif, percent_overlap)
            datum = [motif, motif2, percent_overlap, 'paratope']
            data.append(datum)
    nodes_epitope = [item for item in penodesdf.id.tolist() if '*' in item][:]
    for motif in nodes_epitope:
        mdf = pedf[pedf.target == motif]
        partners = mdf.source
        for motif2 in nodes_epitope:
            mdf2 = pedf[pedf.target == motif2]
            partners2 = mdf2.source
            intersect = set(partners) & set(partners2)
            percent_overlap = round(len(intersect)/float(len(partners))*100, 1)
            # print(motif, percent_overlap)
            datum = [motif, motif2, percent_overlap, 'epitope']
            data.append(datum)
    rnodes_paratope = [item for item in rnodesdf.id.tolist() if '*' not in item][:]
    for motif in rnodes_paratope:
        mdf = rdf[rdf.source == motif]
        partners = mdf.target
        for motif2 in nodes_paratope:
            mdf2 = rdf[rdf.source == motif2]
            partners2 = mdf2.target
            intersect = set(partners) & set(partners2)
            percent_overlap = round(len(intersect)/float(len(partners))*100, 1)
            # print(motif, percent_overlap)
            datum = [motif, motif2, percent_overlap, 'rparatope']
            data.append(datum)
    rnodes_epitope = [item for item in rnodesdf.id.tolist() if '*' in item][:]
    for motif in rnodes_epitope:
        mdf = rdf[rdf.target == motif]
        partners = mdf.source
        for motif2 in nodes_epitope:
            mdf2 = rdf[rdf.target == motif2]
            partners2 = mdf2.source
            intersect = set(partners) & set(partners2)
            percent_overlap = round(len(intersect)/float(len(partners))*100, 1)
            # print(motif, percent_overlap)
            datum = [motif, motif2, percent_overlap, 'repitope']
            data.append(datum)
    colnames = ['motif1','motif2', 'percent_overlap', 'motif_source']
    print(colnames)
    outdf = pd.DataFrame(data, columns=colnames)
    outname = peinfile.split('.')[0] + '_paratope_cross_density.csv'
    print(outname)
    outdf.to_csv(outname, index=False)




# run stuff
# cross_reactivity_paratope()
# cross_reactivity_paratope_topn()
cross_reactivity_density()
