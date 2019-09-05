# import stuf
import sys
import os
from find_files import find_files as fifi
import pandas as pd
from mpi4py import MPI
import numpy as np



def cross_reactivity_density_paratope_epitope_ppi():
    '''
    prep neat data for cross reactivity density plots both paratope and epitope
    :return:
    '''
    # edge_files = fifi('abdb_outfiles_2019', 'internet_edges.csv')
    # node_files = fifi('abdb_outfiles_2019', 'internet_nodes.csv')
    # infiles = edge_files  + node_files
    # print(infiles)
    # sys.exit()
    peinfile = 'abdb_outfiles_2019/ppi_internet_edges.csv'
    rinfile = 'abdb_outfiles_2019/downsampled_ppi_internet_edges.csv'
    peinfile_nodes = 'abdb_outfiles_2019/ppi_internet_nodes.csv'
    rinfile_nodes = 'abdb_outfiles_2019/downsampled_ppi_internet_nodes.csv'
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
            datum = [motif, motif2, percent_overlap, 'ppimotif']
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
            datum = [motif, motif2, percent_overlap, 'ppimotifpartner']
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
            datum = [motif, motif2, percent_overlap, 'dsppimotif']
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
            datum = [motif, motif2, percent_overlap, 'dsppimotifpartner']
            data.append(datum)
    colnames = ['motif1','motif2', 'percent_overlap', 'motif_source']
    print(colnames)
    outdf = pd.DataFrame(data, columns=colnames)
    outname = peinfile.split('.')[0] + '_paratope_cross_density_ppi.csv'
    print(outname)
    outdf.to_csv(outname, index=False)



def cross_reactivity_density_paratope_epitope_ppi_mpi(nodefile, edgefile, sourcetag, targettag):
    '''
    prep neat data for cross reactivity density plots both paratope and epitope
    retrofit for ppi usage
    :return:
    '''
    # edge_files = fifi('abdb_outfiles_2019', 'internet_edges.csv')
    # node_files = fifi('abdb_outfiles_2019', 'internet_nodes.csv')
    # infiles = edge_files  + node_files
    # print(infiles)
    # sys.exit()
    peinfile = 'abdb_outfiles_2019/ppi_internet_edges.csv'
    rinfile = 'abdb_outfiles_2019/downsampled_ppi_internet_edges.csv'
    peinfile_nodes = 'abdb_outfiles_2019/ppi_internet_nodes.csv'
    rinfile_nodes = 'abdb_outfiles_2019/downsampled_ppi_internet_nodes.csv'
    penodesdf = pd.read_csv(peinfile_nodes)
    rnodesdf = pd.read_csv(rinfile_nodes)
    pedf = pd.read_csv(peinfile).iloc[:]
    rdf = pd.read_csv(rinfile)
    print(pedf.head())
    data = []
    data2 = []
    n = 4
    nodes_paratope = [item for item in penodesdf.id.tolist() if '*' not in item][:]
    nodes_epitope = [item for item in penodesdf.id.tolist() if '*' in item][:]
    # chunk the list
    chunks = np.array_split(nodes_paratope, n)
    chunks2 = np.array_split(nodes_epitope, n)
    print(len(chunks[0]), len(nodes_paratope))
    # scatter the params
    comm = MPI.COMM_WORLD
    print(comm.Get_size())
    if comm.rank == 0:
        params = chunks
        params2 = chunks2
    else:
        params = None
        params2 = None
    params = comm.scatter(params, root=0)
    params2 = comm.scatter(params2, root=0)
    outdir  = 'supfig12outs'
    # clear outdir before making a new one
    os.rmdir(outdir)
    os.mkdir(outdir)
    outname = outdir + '/' + nodefile.split('.')[0] + 'rep%s_%s' % (comm.rank, sourcetag) + '.csv'
    outname2 = outdir + '/' + nodefile.split('.')[0] + 'rep%s_%s' % (comm.rank, targettag) + '.csv'
    print(outname)
    print(outname2)
    print(params, comm.rank)
    print(params2, comm.rank)
    for motif in params:
        mdf = pedf[pedf.source == motif]
        partners = mdf.target
        for motif2 in nodes_paratope:
            mdf2 = pedf[pedf.source == motif2]
            partners2 = mdf2.target
            intersect = set(partners) & set(partners2)
            percent_overlap = round(len(intersect)/float(len(partners))*100, 1)
            # print(motif, percent_overlap)
            datum = [motif, motif2, percent_overlap, 'ppimotif']
            data.append(datum)
    colnames = ['motif1','motif2', 'percent_overlap', 'motif_source']
    outdf1 = pd.DataFrame(data, columns=colnames)
    outdf1.to_csv(outname, index=False)
    for motif in params2:
        mdf = pedf[pedf.target == motif]
        partners = mdf.source
        for motif2 in nodes_epitope:
            mdf2 = pedf[pedf.target == motif2]
            partners2 = mdf2.source
            intersect = set(partners) & set(partners2)
            percent_overlap = round(len(intersect)/float(len(partners))*100, 1)
            # print(motif, percent_overlap)
            datum = [motif, motif2, percent_overlap, 'ppimotifpartner']
            data.append(datum)
    outdf2 = pd.DataFrame(data2, columns=colnames)
    outdf2.to_csv(outname2, index=False)



# run stuff
# cross_reactivity_density_paratope_epitope_ppi()
cross_reactivity_density_paratope_epitope_ppi_mpi('a', 'b', 'c', 'd')