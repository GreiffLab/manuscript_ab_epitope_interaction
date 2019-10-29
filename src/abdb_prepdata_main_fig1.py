# data preparation workflow
# use pdb files in NR_LH_Protein_Martin
# identify intearcting residues by a cutoff: 5 (or any other number)
# output a csv file.

# import stuff
from abdb import *
import sys
import os
from find_files import find_files
import numpy as np

# create outdir
outpath = 'abdb_outfiles_2019'
if os.path.isfile(outpath) == False:
    os.system('mkdir %s' % outpath)
#define a cutoff
cutoff = 5


# examine median resolution
def get_median_resolution():
    '''
    get median resolution in the final dataset
    :return:
    '''
    infile = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_angle_bnaber_phil_pc.csv'
    df = pd.read_csv(infile)
    print(df.head())
    pdbids = df.pdbid.unique()
    print(len(pdbids))
    resolutions =[]
    for pdbid in pdbids:
        pdbfile = '/Users/rahmadakbar/greifflab/aims/aimugen/datasets/NR_LH_Protein_Martin/' + pdbid + '.pdb'
        contents = open(pdbfile).read().splitlines()
        for content in contents[:10]:
            if 'RESOLUTION' in content:
                parts = content.split()
                resolution = float(parts[-1])
                resolutions.append(resolution)
    mean_resolution = round(sum(resolutions)/len(resolutions),2)
    median_resolution = np.median(resolutions)
    print('Median resolution %s, total structures %s' % (median_resolution, len(resolutions)))


# #start
# # get pdb with single antigen
# single_antigens = get_single_antigens()
# # sort by antibody
# absorted, agsorted = get_residue_pairs_ab2(single_antigens[:], outpath, cutoff)
# # account for inserted residues
# abinsert = get_unique_abresnumi(absorted,outpath)
# # add segments based on Martin numbering
# absegment = add_segments(abinsert, outpath)
# # add shift
# abshift = add_abshift(absegment,outpath)
# # add shifft loop wise
# abshiftloop = add_abshiftl(abshift,outpath)
# #get gap patterns data
# gap_patterns = get_numgaps_segment(abshiftloop,outpath)
# # add shift to ag
# add_agshift(abshiftloop, outpath)

# ## per segment run
# segment_files = find_files(outpath, 'segment.csv')
# segment_files = [item for item in segment_files if 'abshift' in item and str(cutoff) in item] # filter for preprocessed
# # make gap dataset
# for segment_file in segment_files[:]:
#     make_gap_dataset(segment_file, outpath)
# batch_add_notationx(segment_files)
# # ouput separate count data for gaps
# notationx_files = find_files('abdb_outfiles_2019', 'notationx.csv')
# notationx_files = [item for item in notationx_files if str(cutoff) in item]
# batch_add_gap_count_data(notationx_files)

get_median_resolution()
