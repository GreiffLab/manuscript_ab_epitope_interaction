
# import stuff
import sys
sys.path.insert(0, 'c:\\anaconda3\\lib\\site-packages')
#import pandas as pd
#import jellyfish

import mpi4py
#from mpi4py import MPI
#import numpy as np


import os
import os.path
from os import path
os.environ['PATH'] += 'PATH=%PATH%;C:\\MyPrograms\\MPI\\bin'
#os.system('python C:/Users/pprobert/Desktop/Main/CurrentZapotect/Zapotec/pdb-tools/pdbtools/pdb_selchain.py')

def get_ppi_interacting_segment():
    '''
    get ppi interacting segments from pdb files in
    :return:
    '''
    pdbdir = 'C:/Users/pprobert/Desktop/TempWorld/PPIoriginalPDBs/'
    infile = 'C:/Users/pprobert/Desktop/TempWorld/interacting_residues/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_resnum Original.csv'
    df = pd.read_csv(infile)
    pdbchainpairs = df.pdbchainpair1.tolist()[:]
    # pdbchainpairs = [item for item in pdbchainpairs if item == '1VXQ_NE']
    total_size = len(pdbchainpairs)
    chunks = np.array_split(pdbchainpairs,4)
    data = pdbchainpairs;
    ## scatter with mpi
    # #comm = MPI.COMM_WORLD
    # size = comm.Get_size()
    # rank = comm.Get_rank()
    #
    # if rank == 0:
    #     data = chunks
    # else:
    #     data = None
    #
    # data = comm.scatter(data, root=0)
    # total_size = len(data)
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
            resnum1 = [int(item) for item in cdf.resnum2.iloc[0].split('-')]
            resnum2 = [int(item) for item in cdf.resnum1.iloc[0].split('-')]
            min1, max1 =  min(resnum1)-1, max(resnum1)+1
            min2, max2 =  min(resnum2)-1, max(resnum2)+1
            pdbfile = pdbdir + '/' + pdbid + '.pdb'
            pdbchainfile1 = pdbdir + '/' + pdbid + '_' + chain1 + '.pdb'
            pdbchainfile2 = pdbdir + '/' + pdbid + '_' + chain2 + '.pdb'
            command1 = 'python C:/Users/pprobert/Desktop/Main/CurrentZapotect/Zapotec/pdb-tools/pdbtools/pdb_selchain.py -%s %s > %s' % (chain1, pdbfile, pdbchainfile1)
            command2 = 'python C:/Users/pprobert/Desktop/Main/CurrentZapotect/Zapotec/pdb-tools/pdbtools/pdb_selchain.py -%s %s > %s' % (chain2, pdbfile, pdbchainfile2)
            os.system(command1)
            os.system(command2)
            segmentdir = 'C:/Users/pprobert/Desktop/TempWorld/NewPPIcutPDBs'
            sfile1 = segmentdir + '/' + pdbid + '_' + chain1 + '_X1_%s%s_' % (chain1, chain2) + '_interacting_segment_%s.pdb' % tag
            sfile2 = segmentdir + '/' + pdbid + '_' + chain2 + '_X2_%s%s_' % (chain1, chain2) + '_interacting_segment_%s.pdb' % tag
            if(not path.exists(sfile1)):
                scommand1 = 'python C:/Users/pprobert/Desktop/Main/CurrentZapotect/Zapotec/pdb-tools/pdbtools/pdb_selres.py -%s:%s %s > %s' % (min1, max1, pdbchainfile1, sfile1)
                scommand2 = 'python C:/Users/pprobert/Desktop/Main/CurrentZapotect/Zapotec/pdb-tools/pdbtools/pdb_selres.py -%s:%s %s > %s' % (min2, max2, pdbchainfile2, sfile2)
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
                #print('str %s done of %s, rank %s' %  (i, total_size, rank))
                # print('str %s of %s' % (i, total_size))
                # print(cdf)
                # sys.exit()
        except:
            except_counter += 1
            print(except_counter)
            print('excepting...')


#get_ppi_interacting_segment();
#sys.exit();



def get_ppi_interacting_segmentMPI():
    '''
    get ppi interacting segments from pdb files in
    :return:
    '''
    pdbdir = 'C:/Users/pprobert/Desktop/TempWorld/PPIoriginalPDBs/'
    infile = 'C:/Users/pprobert/Desktop/TempWorld/interacting_residues/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_resnum Original.csv'
    df = pd.read_csv(infile)
    pdbchainpairs = df.pdbchainpair1.tolist()[:]
    # pdbchainpairs = [item for item in pdbchainpairs if item == '1VXQ_NE']
    total_size = len(pdbchainpairs)
    chunks = np.array_split(pdbchainpairs,4)
    #data = pdbchainpairs;
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
    # end MPI
    # for i,pdbchainpair in enumerate(pdbchainpairs):
    except_counter = 0
    for i,pdbchainpair in enumerate(data):
        #try:
            parts  = pdbchainpair.split('_')
            pdbid = parts[0]
            chain1 = parts[1][0]
            chain2 =  parts[1][1]
            tag = 'inter'
            if chain1 == chain2:
                tag = 'intra'
            cdf = df[df.pdbchainpair1 == pdbchainpair]
            resnum1 = [int(item) for item in cdf.resnum2.iloc[0].split('-')]
            resnum2 = [int(item) for item in cdf.resnum1.iloc[0].split('-')]
            min1, max1 =  min(resnum1)-1, max(resnum1)+1
            min2, max2 =  min(resnum2)-1, max(resnum2)+1
            pdbfile = pdbdir + '/' + pdbid + '.pdb'
            pdbchainfile1 = pdbdir + '/' + pdbid + '_' + chain1 + '.pdb'
            pdbchainfile2 = pdbdir + '/' + pdbid + '_' + chain2 + '.pdb'
            command1 = 'python C:/Users/pprobert/Desktop/Main/CurrentZapotect/Zapotec/pdb-tools/pdbtools/pdb_selchain.py -%s %s > %s' % (chain1, pdbfile, pdbchainfile1)
            command2 = 'python C:/Users/pprobert/Desktop/Main/CurrentZapotect/Zapotec/pdb-tools/pdbtools/pdb_selchain.py -%s %s > %s' % (chain2, pdbfile, pdbchainfile2)
            os.system(command1)
            os.system(command2)
            segmentdir = 'C:/Users/pprobert/Desktop/TempWorld/NewPPIcutPDBs'
            sfile1 = segmentdir + '/' + pdbid + '_' + chain1 + '_X1_%s%s_' % (chain1, chain2) + '_interacting_segment_%s.pdb' % tag
            sfile2 = segmentdir + '/' + pdbid + '_' + chain2 + '_X2_%s%s_' % (chain1, chain2) + '_interacting_segment_%s.pdb' % tag
            scommand1 = 'python C:/Users/pprobert/Desktop/Main/CurrentZapotect/Zapotec/pdb-tools/pdbtools/pdb_selres.py -%s:%s %s > %s' % (min1, max1, pdbchainfile1, sfile1)
            scommand2 = 'python C:/Users/pprobert/Desktop/Main/CurrentZapotect/Zapotec/pdb-tools/pdbtools/pdb_selres.py -%s:%s %s > %s' % (min2, max2, pdbchainfile2, sfile2)
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
        #except:
        #    except_counter += 1
        #    print(except_counter)
        #    print('excepting...')


#get_ppi_interacting_segmentMPI();
#sys.exit();








#I need to tell where are the librairies (non-admin priviledge honour to do it)
sys.path.insert(0, 'c:\\anaconda3\\lib\\site-packages')
import Bio
from Bio import PDB

#----------- 1: takes all the PDB in one file and generates a text file with all ramachandran angles of all residues of all chains --------------
# https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/
def folderRamachandran(folder, outputFileNameTxt):
    import os
    sys.path.insert(0, 'c:\\anaconda3\\lib\\site-packages')
    import Bio
    from Bio import PDB
    fwrite= open(folder + outputFileNameTxt,"w+")
    for filename in os.listdir(folder):
        if filename.endswith(".pdb"):
            print(filename)
            for model in Bio.PDB.PDBParser().get_structure("XXXX", folder + filename):
                for chain in model:
                    polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
                    for poly_index, poly in enumerate(polypeptides):
                        phi_psi = poly.get_phi_psi_list()
                        for res_index, residue in enumerate(poly):
                            res_name = "%s\t%i" % (residue.resname, residue.id[1])
                            #print(filename, res_name, "\t", phi_psi[res_index][0], "\t", phi_psi[res_index][1])
                            fwrite.write(str(filename) + "\t" + res_name + "\t" + str(chain._id) + "\t" + str(phi_psi[res_index][0]) + "\t" + str(phi_psi[res_index][1]) + "\n")
    fwrite.close()

#folder = 'C:/Users/pprobert/Desktop/TempWorld/interacting_residues/';
folder = 'C:/Users/pprobert/Desktop/TempWorld/NewPPIcutPDBs/';
outputFileName = 'AnalyzedAllFolderPPIcut.txt'
folderRamachandran(folder, outputFileName)

#----------- 2: Detects PDB files with multiple models --------------
## opening excel file with motifs and its positions, ane returns the PDB with the positions to fetch.
sys.path.insert(0, 'c:\\anaconda3\\lib\\site-packages')
import Bio
from Bio import PDB
import csv;
folder = 'C:/Users/pprobert/Desktop/TempWorld/PPIoriginalPDBs/';
folderXls = 'C:/Users/pprobert/Desktop/TempWorld/interacting_residues/';
csv_file = open(folderXls + 'threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_resnum.csv')
csv_reader = csv.reader(csv_file, delimiter=',')
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
def testPDBs():
    alreadyDone = {''};
    line_count = 0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)
        for row in csv_reader:
            if line_count == 0:
                line_count += 1;
                #print(f'Column names are {", ".join(row)}')
            else:
                #print(row[1], "\t", row[-2]);  # , "\n")
                PDBname = row[1];
                PDBonly = PDBname[0:4]  #note: when slicing, the second index is excluded
                if(not PDBonly in alreadyDone):
                    alreadyDone.add(PDBonly);
                    #print("reading " + PDBonly + '.pdb');
                    Models = Bio.PDB.PDBParser().get_structure("WhateverID", folder + PDBonly + '.pdb');
                    if(len(Models) > 1):
                        print (PDBname + 'is multiModels');
                    else:
                        print(PDBname + 'ok');

#testPDBs();


#https://stackoverflow.com/questions/12760271/how-do-i-convert-the-three-letter-amino-acid-codes-to-one-letter-code-with-pytho
d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def shorten(x):
    if len(x) % 3 != 0:
        raise ValueError('Input length should be a multiple of three')

    y = ''
    for i in range(int(len(x)/3)):
            y += d[x[3*i:3*i+3]]
    return y

shorten('CYSARG')

import os.path
from os import path
import re

def folderFileAndPos(filename, chainID, positions, checksequence, openedOutputFile, ferrors, expectedSequence2):
    #fwrite=open(outputfilename, "a+")
    foundSequence = "";
    result = "";
    if(not path.exists(folder + filename)):
        ferrors.write("ERROR with file " + folder + filename + "Could not open\n");
    else:
        for model in Bio.PDB.PDBParser().get_structure("WhateverID", folder + filename):
            for chain in model:
                #print("chain", chain._id);
                if (chain._id == chainID):
                    polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
                    for poly_index, poly in enumerate(polypeptides):
                        phi_psi = poly.get_phi_psi_list()
                        for res_index, residue in enumerate(poly):
                            #print(residue.id[1], residue.resname, checksequence, positions); #this is the ID
                            #print(residue.id);
                            #if(residue.id[2]):
                                #print(residue.id[2]);
                            #print(str(residue.id[1]) + residue.id[2]);
                            if(((str(str(residue.id[1]) + residue.id[2])) in positions) or ((str(residue.id[1]) in positions) and (residue.id[2] == ' '))): #Fuck Python str(residue.id[1]) + residue.id[2]) in positions doesnt work for mormal numbers
                                #print("found", (str(residue.id[1]) + residue.id[2]));
                                #if(re.search('[a-zA-Z]', str(residue.id[1]))):
                                #    print("detected insertion ", residue.id[1]);
                                res_name = "%s\t%i" % (residue.resname, residue.id[1])
                                foundSequence += d[res_name[0:3]] #for whatever reason it has more characters
                                #print(filename, chain._id, res_name, "\t", phi_psi[res_index][0], "\t", phi_psi[res_index][1])
                                result += str(filename) + "\t" + chain._id + "\t" + res_name + "\t" + str(phi_psi[res_index][0]) + "\t" + str(phi_psi[res_index][1]) + "\n";
        if(((len(checksequence) > 0) and (''.join(sorted(checksequence)) != ''.join(sorted(foundSequence))))  and ((len(expectedSequence2) > 0) and (''.join(sorted(expectedSequence2)) != ''.join(sorted(foundSequence))))) :
            ferrors.write("ERROR in file " + filename + ", the requested interaction motif has different residues. Requested " + checksequence + " or " +  expectedSequence2 + " while found is " + foundSequence+ "\n")
        else:
            openedOutputFile.write(result);

folder = 'C:/Users/pprobert/Desktop/TempWorld/PPIcutPDBs/';
filename = "5UZ4_C_CJ__interacting_segment_inter.pdb"
chainID = 'C'
positions = ('94','95');
seq = 'AG';
outputFileName = folder + '../AnalyzedMotifsPointsOnly.txt';
fwrite = open(outputFileName, "w+");
ferrors = open(outputFileName + "errors.txt", "w+");
#folderFileAndPos(filename, chainID, positions, seq, fwrite, ferrors);
ferrors.close();
fwrite.close();


def pickOnlyResidues(folder, csv_file, outputfilename, fileErrors):
    badFiles = {'1A1T', '1AIY', '1AQ5', '1CNP', '1EJP', '1EJQ', '1EWW', '1HZ8', '1JM7', '1JUN', '1L1I', '1M7L', '1M8O',
                '1MNT', '1MV4', '1MVZ', '1N4I', '1NTC', '1OEI', '1OHH', '1OLH', '1PET', '1PZR', '1QEY', '1RSO', '1SAE',
                '1SAF', '1SAH', '1SAJ', '1SMZ', '1TF3', '1WJB', '1WJD', '1WJF', '1X9V', '1Y74', '1Y76', '1YRQ', '1ZLL',
                '1ZXA', '2AIY', '2BBI', '2CT1', '2DLQ', '2DRN', '2E34', '2EC7', '2EZY', '2EZZ', '2H8W', '2HAC', '2HYN',
                '2J0Z', '2J10', '2J11', '2JPT', '2JQ7', '2JTT', '2K1N', '2K3F', '2K6S', '2K7O', '2K8M', '2KA1', '2KA2',
                '2KBY', '2KI6', '2KIX', '2KJ1', '2KKG', '2KKO', '2KPE', '2KPF', '2KQT', '2KU1', '2KWX', '2KYV', '2L14',
                '2L48', '2L49', '2L5G', '2L6Y', '2L6Z', '2L9N', '2LJB', '2LJC', '2LMN', '2LMO', '2LMP', '2LMQ', '2LNH',
                '2LQL', '2LY0', '2LZ3', '2LZ4', '2M3B', '2M56', '2M59', '2M89', '2M8R', '2MET', '2MIC', '2MJO', '2MJZ',
                '2MK9', '2MKA', '2ML1', '2MN6', '2MPZ', '2MUV', '2MUW', '2MXU', '2N1F', '2N1T', '2N27', '2N3U', '2N3W',
                '2N70', '2N90', '2N9B', '2NB1', '2OTK', '2Q4G', '2Q50', '2RLF', '3JCK', '3SAK', '4B2S', '4TVX', '4WT8',
                '5AAS', '5FM1', '5HUZ', '5IJN', '5JTN', '5K5G', '5KK3', '5N5A', '5NZT', '5TN2', '5UZ9', '5V4U', '5V7Z',
                '5W0S', '5YAM', '6AMW', '6AZ0', '6BZL', '6DLN'};
    ## opening excel file with motifs and its positions, ane returns the PDB with the positions to fetch.
    import csv;
    csv_reader = csv.reader(csv_file, delimiter=',')
    fwrite = open(outputfilename, "w+");
    fwrite.write("PDB\tchain\tresidue\tposition\tphi\tpsi\n");
    ferrors = open(fileErrors, "w+");
    line_count = 0
    for row in csv_reader:
        if line_count > 100000:
            return
        if line_count == 0:
            line_count += 1;
            #print(f'Column names are {", ".join(row)}')
        else:
            line_count += 1;
            #print(row[1], "\t", row[-2]);  # , "\n")
            PDBname = row[-1];
            PDBonly = PDBname[0:4]  #note: when slicing, the second index is excluded
            if(not PDBonly in badFiles):
                # bla
                chain1 = PDBname[5];
                chain2 = PDBname[6];
                # bla
                toCut1 = row[20];
                expectedSequence1 = row[2];
                type =  row[6][0:5]; #inter or intra
                listPos1 = toCut1.split('-')
                toCut2 = row[21];
                expectedSequence2 = row[12];
                listPos2 = toCut2.split('-')
                if(len(expectedSequence1) != len(listPos1)):
                    if(len(expectedSequence1) != len(listPos2)):
                        print("Inconsistent both ways" + expectedSequence1 + " " + str(len(expectedSequence1)) + " " + str(len(listPos1)) + " " + str(toCut1) + " " + expectedSequence2 + str(len(expectedSequence2)));
                    store = expectedSequence1;
                    expectedSequence1 = expectedSequence2;
                    expectedSequence2 = store;


                #print(expectedSequence1, listPos1);
                #print(expectedSequence2, listPos2);
                filename1 = PDBonly + "_" + chain1 + "_X2_" + chain2 + chain1 + "__interacting_segment_" + type + ".pdb";
                folderFileAndPos(filename1, chain1, listPos1, expectedSequence1, fwrite, ferrors, expectedSequence2);
                # bla
                filename2 = PDBonly + "_" + chain2 + "_X1_" + chain2 + chain1 + "__interacting_segment_" + type + ".pdb";
                folderFileAndPos(filename2, chain2, listPos2, expectedSequence2, fwrite, ferrors, expectedSequence1);
    fwrite.close();
    ferrors.close();


folder = 'C:/Users/pprobert/Desktop/TempWorld/NewPPIcutPDBs/';
folderXls = 'C:/Users/pprobert/Desktop/TempWorld/interacting_residues/';
csv_file = open(folderXls + 'threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_resnum_flipped.csv')
outputFileName = folder + '../AnalyzedMotifsPointsOnly.txt'
fileErrors = folder + '../AnalyzedMotifsPointsOnlyERRORS.txt'
pickOnlyResidues(folder, csv_file, outputFileName, fileErrors);


def folderFileAndPosSimpler(filename, chainID, positions, checksequence, openedOutputFile, ferrors):
    #fwrite=open(outputfilename, "a+")
    foundSequence = "";
    result = "";
    if(not path.exists(folder + filename)):
        ferrors.write("ERROR with file " + folder + filename + "Could not open\n");
    else:
        for model in Bio.PDB.PDBParser().get_structure("WhateverID", folder + filename):
            for chain in model:
                #print("chain", chain._id);
                if (chain._id == chainID):
                    polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
                    for poly_index, poly in enumerate(polypeptides):
                        phi_psi = poly.get_phi_psi_list()
                        for res_index, residue in enumerate(poly):
                            #print(residue.id[1], residue.resname, checksequence, positions); #this is the ID
                            #print(residue.id);
                            #if(residue.id[2]):
                                #print(residue.id[2]);
                            #print(str(residue.id[1]) + residue.id[2]);
                            if(((str(str(residue.id[1]) + residue.id[2])) in positions) or ((str(residue.id[1]) in positions) and (residue.id[2] == ' '))): #Fuck Python str(residue.id[1]) + residue.id[2]) in positions doesnt work for mormal numbers
                                #print("found", (str(residue.id[1]) + residue.id[2]));
                                #if(re.search('[a-zA-Z]', str(residue.id[1]))):
                                #    print("detected insertion ", residue.id[1]);
                                res_name = "%s\t%i" % (residue.resname, residue.id[1])
                                foundSequence += d[res_name[0:3]] #for whatever reason it has more characters
                                #print(filename, chain._id, res_name, "\t", phi_psi[res_index][0], "\t", phi_psi[res_index][1])
                                result += str(filename) + "\t" + chain._id + "\t" + res_name + "\t" + str(phi_psi[res_index][0]) + "\t" + str(phi_psi[res_index][1]) + "\n";
        if((len(checksequence) > 0) and (''.join(sorted(checksequence)) != ''.join(sorted(foundSequence)))) : #cannot get them in the same order. so sort them and compare
            ferrors.write("ERROR in file " + filename + ", the requested interaction motif has different residues. Requested " + checksequence + " while found is " + foundSequence+ "\n")
        else:
            openedOutputFile.write(result);

def ABDBpickOnlyResidues(folder, csv_file, outputfilename, fileErrors):
    badFiles = {''};
    ## opening excel file with motifs and its positions, ane returns the PDB with the positions to fetch.
    import csv;
    csv_reader = csv.reader(csv_file, delimiter=',')
    fwrite = open(outputfilename, "w+");
    fwrite.write("PDB\tchain\tresidue\tposition\tphi\tpsi\n");
    ferrors = open(fileErrors, "w+");
    line_count = 0
    for row in csv_reader:
        if line_count > 100000:
            return
        if line_count == 0:
            line_count += 1;
            #print(f'Column names are {", ".join(row)}')
        else:
            line_count += 1;
            #print(row[1], "\t", row[-2]);  # , "\n")
            PDBname = row[0];
            PDBonly = PDBname[0:4]  #note: when slicing, the second index is excluded
            if(not PDBonly in badFiles):
                chain1 = row[1];
                toCut1 = row[7];
                expectedSequence1 = row[3];
                listPos1 = toCut1.split('-')
                #print(expectedSequence1, listPos1);
                filename1 = PDBname + ".pdb";
                folderFileAndPosSimpler(filename1, chain1, listPos1, expectedSequence1, fwrite, ferrors);
    fwrite.close();
    ferrors.close();


folder = 'C:/Users/pprobert/Desktop/TempWorld/ABDBoriginalPDBs/';
folderXls = 'C:/Users/pprobert/Desktop/TempWorld/';
csv_file = open(folderXls + 'respairs_segment_notationx_len_merged_angle_bnaber_phil_pc.csv')
outputFileName = folder + '../AnalyzedABDBMotifsPointsOnly.txt'
fileErrors = folder + '../AnalyzedABDBMotifsPointsOnlyERRORS.txt'
#ABDBpickOnlyResidues(folder, csv_file, outputFileName, fileErrors);





def remove_insertion(positionAsString):
    Myres = "";
    listNumbers = {'0', '1', '2', '3', '4', '5', '6', '7' , '8', '9'};
    for letter in positionAsString:
        if(letter in listNumbers):
            Myres = Myres + letter;
        else:
            return Myres;
    return Myres;

#remove_insertion('123C1')


def ABDBpickFullMotifsResidues(folder, csv_file, outputfilename, fileErrors):
    badFiles = {''};
    ## opening excel file with motifs and its positions, ane returns the PDB with the positions to fetch.
    import csv;
    csv_reader = csv.reader(csv_file, delimiter=',')
    fwrite = open(outputfilename, "w+");
    fwrite.write("PDB\tchain\tresidue\tposition\tphi\tpsi\n");
    ferrors = open(fileErrors, "w+");
    line_count = 0
    for row in csv_reader:
        #if(line_count > 20):
        #    return;
        #if line_count > 100000:
        #    return
        if line_count == 0:
            line_count += 1;
            #print(f'Column names are {", ".join(row)}')
        else:
            line_count += 1;
            #print(row[1], "\t", row[-2]);  # , "\n")
            PDBname = row[0];
            PDBonly = PDBname[0:4]  #note: when slicing, the second index is excluded
            if(not PDBonly in badFiles):
                chain1 = row[1];
                toCut1 = row[7];
                expectedSequence1 = "";
                listPos1 = toCut1.split('-')
                #print(listPos1);
                Mmin = +1e6;
                Mmax = -100;
                for point in listPos1:
                    Mmin = min(Mmin, float(remove_insertion(point)));
                    Mmax = max(Mmax, float(remove_insertion(point)));
                newListPos = list(range(int(max(0,Mmin-1)), int(Mmax+2))); #in range, the last value is nto tken
                #newListPosInString = map(str, newListPos); didnt work
                newListPosInString = [str(i) for i in newListPos];
                #print(newListPosInString);
                #print(listPos1 + "\n => " + newListPosInString + "\n");
                #print(expectedSequence1, listPos1);
                filename1 = PDBname + ".pdb";
                folderFileAndPosSimpler(filename1, chain1, newListPosInString, expectedSequence1, fwrite, ferrors);
    fwrite.close();
    ferrors.close();



folder = 'C:/Users/pprobert/Desktop/TempWorld/ABDBoriginalPDBs/';
folderXls = 'C:/Users/pprobert/Desktop/TempWorld/';
csv_file = open(folderXls + 'respairs_segment_notationx_len_merged_angle_bnaber_phil_pc.csv')
outputFileName = folder + '../AnalyzedABDBALLMotifsOnly.txt'
fileErrors = folder + '../AnalyzedABDBALLMotifsOnlyERRORS.txt'
ABDBpickFullMotifsResidues(folder, csv_file, outputFileName, fileErrors);


#
#
#
#
#
# ## now opens the PDBs
# folder = 'C:/Users/pprobert/Desktop/TempWorld/interacting_residues/';
# fwrite= open(folder + 'Analyzed.txt',"w+");
#
# csv_file = open('C:/Users/pprobert/Desktop/TempWorld/interacting_residues/threedid_no_iglike_notationx_merged_maxgap7_maxlen300_paired_resnum.csv')
# csv_reader = csv.reader(csv_file, delimiter=',')
# line_count = 0
# for row in csv_reader:
#     if line_count == 0:
#         print(f'Column names are {", ".join(row)}')
#     else:
#         print(row[1], "\t", row[-2]);  # , "\n")
#         toCut = row[-1];
#         listPos = toCut.split('-')
#         for pos in listPos:
#             print(pos, " ");
#     line_count += 1
#
#     import os
#
#     fwrite = open("C:/Users/pprobert/Desktop/Main/pythonProjects/tempWorld/Analyzed.txt", "w+")
#     for filename in os.listdir("C:/Users/pprobert/Desktop/TempWorld/interacting_residues/interacting_residues"):
#         if filename.endswith(".pdb"):
#             # print(filename)
#             for model in Bio.PDB.PDBParser().get_structure("5FGB_CDR-L1",
#                                                            "C:/Users/pprobert/Desktop/TempWorld/interacting_residues/interacting_residues/" + filename):
#                 for chain in model:
#                     polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
#                     for poly_index, poly in enumerate(polypeptides):
#                         phi_psi = poly.get_phi_psi_list()
#                         for res_index, residue in enumerate(poly):
#                             res_name = "%s\t%i" % (residue.resname, residue.id[1])
#                             # print(filename, res_name, "\t", phi_psi[res_index][0], "\t", phi_psi[res_index][1])
#                             fwrite.write(
#                                 str(filename) + "\t" + res_name + "\t" + str(phi_psi[res_index][0]) + "\t" + str(
#                                     phi_psi[res_index][1]) + "\n")
#     fwrite.close()
#
#
#
#
# def identifyDoubleChains(folder):
#     import os
#     sys.path.insert(0, 'c:\\anaconda3\\lib\\site-packages')
#     import Bio
#     from Bio import PDB
#     for filename in os.listdir(folder):
#         if filename.endswith(".pdb"):
#             print(filename)
#             for model in Bio.PDB.PDBParser().get_structure("XXXX", folder + filename):
#                 for chain in model:
#                     lastIndex = -1;
#                     polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
#                     for poly_index, poly in enumerate(polypeptides):
#                         phi_psi = poly.get_phi_psi_list()
#                         for res_index, residue in enumerate(poly):
#                             if (res_index in positions):
#                                 if (res_index) < lastIndex:
#                                     print(filename)
#                                     return
#                             res_name = "%s\t%i" % (residue.resname, residue.id[1])
#                             #print(filename, res_name, "\t", phi_psi[res_index][0], "\t", phi_psi[res_index][1])
#                             #fwrite.write(str(filename) + "\t" + res_name + "\t" + str(chain._id) + "\t" + str(phi_psi[res_index][0]) + "\t" + str(phi_psi[res_index][1]) + "\n")
#     #fwrite.close()

#----------- 1: takes all the PDB in one file and generates a text file with all ramachandran angles of all residues of all chains --------------
# https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/
def folderRamachandranPPI(folder, outputFileNameTxt):
    badFiles = {'1A1T', '1AIY', '1AQ5', '1CNP', '1EJP', '1EJQ', '1EWW', '1HZ8', '1JM7', '1JUN', '1L1I', '1M7L', '1M8O',
                '1MNT', '1MV4', '1MVZ', '1N4I', '1NTC', '1OEI', '1OHH', '1OLH', '1PET', '1PZR', '1QEY', '1RSO', '1SAE',
                '1SAF', '1SAH', '1SAJ', '1SMZ', '1TF3', '1WJB', '1WJD', '1WJF', '1X9V', '1Y74', '1Y76', '1YRQ', '1ZLL',
                '1ZXA', '2AIY', '2BBI', '2CT1', '2DLQ', '2DRN', '2E34', '2EC7', '2EZY', '2EZZ', '2H8W', '2HAC', '2HYN',
                '2J0Z', '2J10', '2J11', '2JPT', '2JQ7', '2JTT', '2K1N', '2K3F', '2K6S', '2K7O', '2K8M', '2KA1', '2KA2',
                '2KBY', '2KI6', '2KIX', '2KJ1', '2KKG', '2KKO', '2KPE', '2KPF', '2KQT', '2KU1', '2KWX', '2KYV', '2L14',
                '2L48', '2L49', '2L5G', '2L6Y', '2L6Z', '2L9N', '2LJB', '2LJC', '2LMN', '2LMO', '2LMP', '2LMQ', '2LNH',
                '2LQL', '2LY0', '2LZ3', '2LZ4', '2M3B', '2M56', '2M59', '2M89', '2M8R', '2MET', '2MIC', '2MJO', '2MJZ',
                '2MK9', '2MKA', '2ML1', '2MN6', '2MPZ', '2MUV', '2MUW', '2MXU', '2N1F', '2N1T', '2N27', '2N3U', '2N3W',
                '2N70', '2N90', '2N9B', '2NB1', '2OTK', '2Q4G', '2Q50', '2RLF', '3JCK', '3SAK', '4B2S', '4TVX', '4WT8',
                '5AAS', '5FM1', '5HUZ', '5IJN', '5JTN', '5K5G', '5KK3', '5N5A', '5NZT', '5TN2', '5UZ9', '5V4U', '5V7Z',
                '5W0S', '5YAM', '6AMW', '6AZ0', '6BZL', '6DLN'};
    import os
    sys.path.insert(0, 'c:\\anaconda3\\lib\\site-packages')
    import Bio
    from Bio import PDB
    fwrite= open(folder + outputFileNameTxt,"w+")
    for filename in os.listdir(folder):
        if filename.endswith(".pdb"):
            print(filename)
            for model in Bio.PDB.PDBParser().get_structure("XXXX", folder + filename):
                for chain in model:
                    polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
                    for poly_index, poly in enumerate(polypeptides):
                        phi_psi = poly.get_phi_psi_list()
                        for res_index, residue in enumerate(poly):
                            res_name = "%s\t%i" % (residue.resname, residue.id[1])
                            #print(filename, res_name, "\t", phi_psi[res_index][0], "\t", phi_psi[res_index][1])
                            fwrite.write(str(filename) + "\t" + res_name + "\t" + str(chain._id) + "\t" + str(phi_psi[res_index][0]) + "\t" + str(phi_psi[res_index][1]) + "\n")
    fwrite.close()

#folder = 'C:/Users/pprobert/Desktop/TempWorld/interacting_residues/';
folder = 'C:/Users/pprobert/Desktop/TempWorld/PPIcutPDBs/';
outputFileName = 'AnalyzedAllFolder.txt'
#folderRamachandran(folder, outputFileName)