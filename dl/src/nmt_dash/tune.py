# param tuning
# # training, split this into a diff code.
from nmt import *
import pandas as pd
import os

infile = '../../dataset/paraepi.tsv'
# infile = '../../dataset/motif_epiparadash.tsv'
input_tensor, target_tensor, inp_lang, targ_lang, max_length_inp, max_length_targ = load_dataset(infile)
print(inp_lang.idx2word)
print(targ_lang.idx2word)
input_tensor_train, input_tensor_val, target_tensor_train, target_tensor_val = train_test_split(input_tensor, target_tensor, test_size=0.2)
inchar_valchar = ''

#make dirs for training checkpoints, losses and test files
if 'motif' in infile:
    sep = ' '
    checkpoint_dir_prefix = './motif_training_checkpoints'
    test_dir = 'motif_test_files'
    loss_dir = 'motif_training_losses'
else:
    sep = ''
    checkpoint_dir_prefix = './training_checkpoints'
    test_dir = 'test_files'
    loss_dir = './training_losses'
tags = ['<start>', '<end>', '<pad>']

for inval, tarval in zip(input_tensor_val,target_tensor_val):
    # print(inval, tarval)
    inchar  = sep.join([inp_lang.idx2word[item] for item in inval if inp_lang.idx2word[item] not in tags])

    valchar  = sep.join([targ_lang.idx2word[item] for item in tarval if targ_lang.idx2word[item] not in tags])
    # print(inchar,'brek', valchar)
    inchar_valchar += '\t'.join([inchar, valchar]) + '\n'

os.system('mkdir %s' %test_dir)
test_outpath = test_dir + '/' + infile.split('/')[-1]
test_outfile = open(test_outpath, 'w')
test_outfile.write(inchar_valchar)
print(test_outfile)
# check length
print(len(input_tensor_train), len(target_tensor_train), len(input_tensor_val), len(target_tensor_val))

#tune embedding dims and unit sizes
epoch_loss = []
embedding_dims = [4**i for i in range(1, 2)]
unit_sizes  = [4**i for i in range(1, 2)]
print(embedding_dims)
for embdim in embedding_dims:
    embedding_dim = embdim
    for usizes in unit_sizes:
        units = usizes
    ### params
        BUFFER_SIZE = len(input_tensor_train)
        BATCH_SIZE = 32
        N_BATCH = BUFFER_SIZE//BATCH_SIZE
        # embedding_dim = 4
        # units = 4
        vocab_inp_size = len(inp_lang.word2idx)
        vocab_tar_size = len(targ_lang.word2idx)

        dataset = tf.data.Dataset.from_tensor_slices((input_tensor_train, target_tensor_train)).shuffle(BUFFER_SIZE)
        dataset = dataset.batch(BATCH_SIZE, drop_remainder=True)
        encoder = Encoder(vocab_inp_size, embedding_dim, units, BATCH_SIZE)
        decoder = Decoder(vocab_tar_size, embedding_dim, units, BATCH_SIZE)
        optimizer = tf.train.AdamOptimizer()

        #check point outfiles
        checkpoint_dir = checkpoint_dir_prefix + '_embdim%s_unit%s' % (embedding_dim, units)
        checkpoint_prefix = os.path.join(checkpoint_dir, "ckpt")
        checkpoint = tf.train.Checkpoint(optimizer=optimizer,
                                         encoder=encoder,
                                         decoder=decoder)

        EPOCHS = 2
        for epoch in range(EPOCHS):
            start = time.time()

            hidden = encoder.initialize_hidden_state()
            total_loss = 0

            for (batch, (inp, targ)) in enumerate(dataset):
                loss = 0
                # print(batch, (inp,targ))
                # sys.exit()

                with tf.GradientTape() as tape:
                    enc_output, enc_hidden = encoder(inp, hidden)

                    dec_hidden = enc_hidden

                    dec_input = tf.expand_dims([targ_lang.word2idx['<start>']] * BATCH_SIZE, 1)

                    # Teacher forcing - feeding the target as the next input
                    for t in range(1, targ.shape[1]):
                        # passing enc_output to the decoder
                        predictions, dec_hidden, _ = decoder(dec_input, dec_hidden, enc_output)

                        loss += loss_function(targ[:, t], predictions)

                        # using teacher forcing
                        dec_input = tf.expand_dims(targ[:, t], 1)

                batch_loss = (loss / int(targ.shape[1]))

                total_loss += batch_loss

                variables = encoder.variables + decoder.variables

                gradients = tape.gradient(loss, variables)

                optimizer.apply_gradients(zip(gradients, variables))

                if batch % 100 == 0:
                    print('Epoch {} Batch {} Loss {:.4f}'.format(epoch + 1,
                                                                 batch,
                                                                 batch_loss.numpy()))
            # saving (checkpoint) the model every 2 epochs
            if (epoch + 1) % 2 == 0:
                checkpoint.save(file_prefix=checkpoint_prefix)

            print('Epoch {} Loss {:.4f}'.format(epoch + 1,
                                                total_loss / N_BATCH))
            print('Time taken for 1 epoch {} sec'.format(time.time() - start))
            print('running embedding dim: %s, units: %s\n' %(embedding_dim, units))
            total_batch_loss = total_loss.numpy()/N_BATCH
            epoch_loss_datum  = [epoch+1, total_batch_loss, time.time()-start, embedding_dim, units]
            epoch_loss.append(epoch_loss_datum)

colnames = ['batch', 'loss', 'time', 'embedding_dim', 'cell_units']
losdf = pd.DataFrame(epoch_loss, columns=colnames)
os.system('mkdir %s' % loss_dir)
out_loss = loss_dir + '/training_loss.csv'
losdf.to_csv(out_loss, index=False)
print(losdf)

