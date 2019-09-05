# import stuff
import sys
import os
import pandas as pd
from abdb_prepdata_main_fig4 import get_context_pattern_top3 as gcpt

# set to display full table
pd.set_option('display.max_column', None)

def split_human_mouse(infile):
    '''
    splits the data to human and mouse for seq dependency net.
    :return:
    '''
    # infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx.csv'
    ref_file = 'abdb_outfiles_2019/respairs_segment_notationx_len_merged_species2.csv'
    df = pd.read_csv(infile)
    refdf = pd.read_csv(ref_file)
    homodf = refdf[refdf.hspecies == 'HOMO SAPIENS']
    musdf = refdf[refdf.hspecies == 'MUS MUSCULUS']
    homo_outdf = df[df.pdbid.isin(homodf.pdbid)]
    mus_outdf = df[df.pdbid.isin(musdf.pdbid)]
    print(homo_outdf.shape)
    print(mus_outdf.shape)
    outname_homo = infile.split('.')[0] + '_homo.csv'
    outname_mus = infile.split('.')[0] + '_mus.csv'
    print(outname_homo)
    print(outname_mus)
    homo_outdf.to_csv(outname_homo, index=False)
    mus_outdf.to_csv(outname_mus, index=False)
    # get edge next
    gcpt(outname_homo)
    gcpt(outname_mus)

## run stuff
split_human_mouse(
    infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_paratope_segment_notationx.csv'

)

split_human_mouse(
    infile = 'abdb_outfiles_2019/respairs_absort_cutoff5_abresnumi_segments_abshift_abshiftl_epitope_segment_notationx'
             '.csv'

)
