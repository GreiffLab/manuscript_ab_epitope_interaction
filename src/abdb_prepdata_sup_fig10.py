# import stuff

import pandas as pd
import sys
import string

# set pd display
pd.set_option('display.max_column', None)

def get_shm_data():
    '''
    check shm from yana, make neat data.
    :return:
    '''
    infile = '../datasets/yana_shm/homo_mus.txt'
    df = pd.read_csv(infile, sep='\t')
    shm_columns = [item for item in df.columns if 'SHMs' in item]
    data = []
    print(shm_columns)
    for i, row in df.iterrows():
        shms = []
        id_dict = {}
        pdbid = row.PDB
        for shm_column in shm_columns:
            idkey = '_'.join(shm_column.split('_')[:2])
            idlen = float(getattr(row, idkey + '_length'))
            chain = idkey[0]
            if getattr(row, shm_column) != '-': # sequences with shm
                contents = getattr(row, shm_column).split(',')
                shmlen = len(contents)
                percent_id = round((idlen-shmlen)/idlen, 2)
                shm_list = [item.split('>')[-1].split(':') + [idkey, percent_id, idlen, shmlen, pdbid, chain] for item
                            in contents]
                if 'V' in shm_column:
                    data += shm_list[:-3]
                elif 'J' in shm_column:
                    data += shm_list[3:]
                else:
                    data += shm_list[:]
                print(shm_list)
            else:
                data += [['-','-', idkey, 1, idlen, 0, pdbid, chain]]
    colnames = ['resnum', 'res', 'chain_gene', 'identity', 'seqlen', 'shmlen', 'pdbid', 'chain']
    outdf = pd.DataFrame(data, columns=colnames)
    tag = infile.split('/')[-1].split('.')[0]
    outname = 'abdb_outfiles_2019/%s_shm_residues.csv' % tag
    outdf.to_csv(outname, index=False)


def get_shm_data_abdb():
    '''
    gets interacting residues from abdb data. incorporates insertion into residue number
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments.csv'
    df = pd.read_csv(infile).iloc[:]
    print(df.head())
    chardict= dict([(item,i+1) for i,item in enumerate(string.ascii_uppercase)])
    chardict[' '] = 0
    abchains = ['L', 'H']
    data = []
    for pdbid in df.pdbid.unique().tolist()[:]:
        pdbfile = '../datasets/NR_LH_Protein_Martin/' + pdbid + '.pdb'
        pdbcontent = open(pdbfile).read().splitlines()
        insertion_dict = {'L':0, 'H':0}
        residues = []
        resnumichains = []
        for line in pdbcontent:
            chain = line[21]
            if line.startswith('ATOM') and chain in abchains:
                resname = line[17:20]
                resnum = int(line[22:26])
                ins = chardict[line[26]]
                resnumichain = resname + str(resnum) + line[26] + chain
                resnumi = resname + str(resnum) + line[26].strip()
                if ins != 0 and resnumichain not in resnumichains:
                    insertion_dict[chain] +=1
                    resnumichains.append(resnumichain)
                insertion = insertion_dict[chain]
                resnum2 = resnum + insertion
                datum = [resname, resnum, ins, resnum2, chain, pdbid, resnumi]
                data.append(datum)

    outname = infile.split('.')[0] + '_adjusted_resnum.csv'
    colnames = ['resname', 'resnum', 'ins', 'resnum2', 'chain', 'pdbid', 'resnumi']
    adjusted_resnums = []
    # df = df.iloc[:100]
    outdf = pd.DataFrame(data=data, columns=colnames)
    for i,row in df.iterrows():
        resnumi = row.abres + str(row.abresnumi)
        pdbid = row.pdbid
        adjdf = outdf[(outdf.pdbid == pdbid) & (outdf.resnumi == resnumi)]
        adjrow = adjdf.iloc[0]
        print('processing: %s/%s' % (i, df.shape[0]))
        adjusted_resnums.append(adjrow.resnum2)
    df['adj_resnum'] = adjusted_resnums
    df.to_csv(outname, index=False)

def process_rescountdf2():
    '''
    add missing residues from 1-131 to rescountdf2
    :return:
    '''
    infile = 'abdb_outfiles_2019/shm_paratope_adjusted_rescountdf2.csv'
    df = pd.read_csv(infile)
    print(df.head())

# run stuff
get_shm_data()
# get_shm_data_abdb()
# process_rescountdf2()

