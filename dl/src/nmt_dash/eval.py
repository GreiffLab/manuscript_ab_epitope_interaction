# test evaluation

from nmt import *
import pandas as pd
import os

# infile = '../../dataset/paraepidash.tsv'
infile = '../../dataset/motif_epiparadash.tsv'
input_tensor, target_tensor, inp_lang, targ_lang, max_length_inp, max_length_targ = load_dataset(infile)
print(input_tensor)
print(inp_lang.idx2word)
print(targ_lang.idx2word)
input_tensor_train, input_tensor_val, target_tensor_train, target_tensor_val = train_test_split(input_tensor, target_tensor, test_size=0.2)

testfile = 'motif_test_files/motif_epiparadash.tsv'
testlines = open(testfile).read().splitlines()
embedding_dims = [4**i for i in range(1, 5)]
unit_sizes  = [4**i for i in range(1, 5)]
evaldata = []
for embdim in embedding_dims:
    embedding_dim = embdim
    for usizes in unit_sizes:
        units = usizes
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
        checkpoint_dir = './motif_training_checkpoints_embdim%s_unit%s' % (embedding_dim, units)
        checkpoint_prefix = os.path.join(checkpoint_dir, "ckpt")
        checkpoint = tf.train.Checkpoint(optimizer=optimizer,
                                         encoder=encoder,
                                         decoder=decoder)

        # # restoring the latest checkpoint in checkpoint_dir
        checkpoint.restore(tf.train.latest_checkpoint(checkpoint_dir))
        for testline in testlines[:]:
            epi, para = testline.split('\t')
            pred_paramotif, epiout =  translate(epi, encoder, decoder, inp_lang, targ_lang, max_length_inp,
                                     max_length_targ)
            tags  = ['<start>', '<end>']
            pred_paramotif = [item for item in pred_paramotif.split() if item not in tags]
            epiout = [item for item in epiout.split() if item not in tags]
            lenpara = len(para.split())
            lenpredpara = len(pred_paramotif)
            lenepi = len(epi)
            evaldatum = [lenpredpara, lenpara,lenepi, embedding_dim, units]
            evaldata.append(evaldatum)

colnames = ['lenpredpara', 'lenpara', 'lenepi', 'embeding_dim', 'units']
evaldf = pd.DataFrame(evaldata, columns=colnames)
print(evaldf)
eval_dir = 'motif_eval_files'
os.system('mkdir %s'%eval_dir)
eval_outpath = eval_dir + '/motif_epipara.csv'
evaldf.to_csv(eval_outpath, index=False)

