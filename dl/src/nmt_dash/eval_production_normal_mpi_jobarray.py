# test evaluation

from nmt import *
import pandas as pd
import os
import jellyfish
from mpi4py import MPI
import sys

# infile = '../../dataset/paraepidash.tsv'
infile = '../../dataset/motif_epiparadash.tsv'
#production dir
prod_dir = '.'
input_tensor, target_tensor, inp_lang, targ_lang, max_length_inp, max_length_targ = load_dataset(infile)
print(input_tensor)
print(inp_lang.idx2word)
print(targ_lang.idx2word)
input_tensor_train, input_tensor_val, target_tensor_train, target_tensor_val = train_test_split(input_tensor, target_tensor, test_size=0.2)

embdims = [1] + [2**i for i in range(1,11)]
lstmunits = [1] + [2**i for i in range(1,11)]
reps = range(1, 11)
rep_emdim_units = [[[(i, j, k) for i in reps] for j in embdims] for k in lstmunits]
rep_emdim_units = sum(rep_emdim_units, [])
rep_emdim_units = sum(rep_emdim_units, [])
rep_emdim_units = rep_emdim_units[:]
embdim_units = [[(i,j) for j in lstmunits] for i in embdims]
embdim_units = sum(embdim_units, [])[:2]

print('total proc: %s' %len(embdim_units))
arrayid = int(sys.argv[1])
repids = range(1,11)

#scatter the params
comm = MPI.COMM_WORLD
print(comm.Get_size())
if comm.rank == 0:
    params = embdim_units
else:
    params = None
params = comm.scatter(params, root=0)
print('rank %s has params %s' % (comm.rank, params))

#for rep, embdim, unit in rep_emdim_units[:]:

for rep in repids[arrayid:arrayid+1]:
    embdim, unit = params
    evaldata = []
    testfile = prod_dir + '/motif_test_files_rep%s/motif_epiparadash.tsv' % rep
    testlines = open(testfile).read().splitlines()
    embedding_dim = embdim
    units = unit
    BUFFER_SIZE = len(input_tensor_train)
    BATCH_SIZE = 32
    N_BATCH = BUFFER_SIZE//BATCH_SIZE
    # embedding_dim = 256
    # units = 256
    vocab_inp_size = len(inp_lang.word2idx)
    vocab_tar_size = len(targ_lang.word2idx)

    encoder = Encoder(vocab_inp_size, embedding_dim, units, BATCH_SIZE)
    decoder = Decoder(vocab_tar_size, embedding_dim, units, BATCH_SIZE)
    optimizer = tf.train.AdamOptimizer()

    #check point outfiles
    checkpoint_dir = '%s/motif_training_checkpoints_rep%s_embdim%s_unit%s' % (prod_dir,rep, embedding_dim, units)
    checkpoint_prefix = os.path.join(checkpoint_dir, "ckpt")
    checkpoint = tf.train.Checkpoint(optimizer=optimizer,
                                     encoder=encoder,
                                     decoder=decoder)

    # # restoring the latest checkpoint in checkpoint_dir
    checkpoint.restore(tf.train.latest_checkpoint(checkpoint_dir))
    for testline in testlines[:10]:
        try:
            epi, para = testline.split('\t')
            pred_paramotif, epiout =  translate(epi, encoder, decoder, inp_lang, targ_lang, max_length_inp,
                                     max_length_targ)
            tags  = ['<start>', '<end>']
            pred_paramotif = [item for item in pred_paramotif.split() if item not in tags]
            str_pred_paramotif = ''.join(pred_paramotif)
            str_para = ''.join(para.split())
            ld = jellyfish.levenshtein_distance(str_pred_paramotif,str_para)
            ldnorm = ld/len(str_pred_paramotif)
            epiout = [item for item in epiout.split() if item not in tags]
            lenpara = len(str_para)
            lenpredpara = len(str_pred_paramotif)
            str_epi = ''.join(epi.split())
            lenepi = len(str_epi)
            evaldatum = [lenpredpara, lenpara,lenepi, rep, embedding_dim, units, ld, ldnorm, str_epi, str_para,
                         str_pred_paramotif]
            evaldata.append(evaldatum)
        except:
            print('skipping this line: %s' % testline)

    colnames = ['lenpredpara', 'lenpara', 'lenepi', 'rep','embeding_dim', 'units', 'ld', 'ldnorm', 'epi', 'para',
                'pre_para']
    evaldf = pd.DataFrame(evaldata, columns=colnames)
    print(evaldf)
    eval_dir = '%s/motif_eval_files' % prod_dir
    os.system('mkdir %s'%eval_dir)
    eval_outpath = eval_dir + '/motif_epipara_rep%s_embdim%s_unit%s.csv' % (rep, embdim, unit)
    evaldf.to_csv(eval_outpath, index=False)

