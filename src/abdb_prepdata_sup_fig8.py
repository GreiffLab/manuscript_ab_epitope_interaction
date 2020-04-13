# import stuff
import pandas as pd
import sys
import os
import numpy as np
import math


def get_angle():
    '''
    gets angels from paratope/epitope
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged.csv'
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
        # print(pdbdf)
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
            try: # some exceptions due to negative residue numbers and insertion in antigens, ignored
                # print(row)
                plen = row.plen
                epitope_len = row.epitope_len
                # print(plen)
                if plen >= 3:
                    # print('hey')
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
                if epitope_len >= 3:
                    paratope = list(row.epitope)
                    abresnumis = row.agresnumiset.split('-')
                    chains = [row.agchain]*len(paratope)
                    centers = []
                    pdbres= ''
                    pseudos = ''
                    for res,chain,abresnumi in zip(paratope, chains, abresnumis):
                        key = ''.join([aadictr[res],chain,abresnumi])
                        # test = [item for item in pdbdict.keys() if item == 'LYSS229']
                        # print(test)
                        # sys.exit()
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
                elif epitope_len <3:
                    e_angle = np.nan
                datum = row.tolist() + [p_angle, e_angle]
                data.append(datum)
            except:
                print('this row failed...ignoring...')
    columns = row.index.tolist() + ['p_angle', 'e_angle']
    outdf = pd.DataFrame(data, columns=columns)
    outname = infile[:-4] + '_angle.csv'
    print(outname)
    outdf.to_csv(outname,index=False)

# run stuff
get_angle()