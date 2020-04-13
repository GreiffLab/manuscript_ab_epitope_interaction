from __future__ import absolute_import, division, print_function

# Import TensorFlow >= 1.10 and enable eager execution
import tensorflow as tf

tf.enable_eager_execution()

import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

import unicodedata
import re
import numpy as np
import os
import time
import pandas as pd
import sys


print(tf.__version__)


def make_tabsep_paraepi():
    '''
    make tab separated para epi file
    :return:
    '''
    infile = '/Users/rahmadakbar/greifflab/aims/aimugen/src/abdb_outfiles/respairs_segment_notationx_merged_angle_len.csv'
    df = pd.read_csv(infile)
    df = df.dropna(subset=['paratope', 'epitope'])
    outcontent = ''
    for i, row in df.iterrows():
        # print(row.epitope, row.paratope)
        content = '\t'.join([row.paratope, row.epitope]) + '\n'
        outcontent += content
    outname = '../../dataset/paraepi.tsv'
    outfile = open(outname, 'w')
    outfile.write(outcontent)
    print(outcontent)

def preprocess_sentence(w):
    '''
    split each character (char is analog to word in a sentence)
    add <start> and <end> token
    :return:
    '''
    return '<start> ' + ' '.join(list(w)) + ' <end>'


def preprocess_sentence_motif(w):
    '''
    split each character (char is analog to word in a sentence)
    add <start> and <end> token
    :return:
    '''
    return '<start> ' + ' '.join(w.split(' ')) + ' <end>'

def create_dataset(infile):
    '''
    split each character (char is analog to word in a sentence)
    add <start> and <end> token
    crates dataset
    :return:
    '''
    lines = open(infile).read().splitlines()
    pairs = []
    for line in lines:
        parts = line.split('\t')
        pair = ['<start> '+' '.join(list(item))+' <end>' for item in parts]
        print(pair, 'seq data')
        pairs.append(pair)
    return pairs


def create_dataset_motif(infile):
    '''
    split each character (char is analog to word in a sentence)
    add <start> and <end> token
    crates dataset
    :return:
    '''
    lines = open(infile).read().splitlines()
    pairs = []
    for line in lines:
        parts = line.split('\t')
        pair = ['<start> '+' '.join(item.split(' '))+' <end>' for item in parts]
        print(pair, 'motif data')
        pairs.append(pair)
    return pairs


class LanguageIndex():
    def __init__(self, lang):
        self.lang = lang
        self.word2idx = {}
        self.idx2word = {}
        self.vocab = set()

        self.create_index()

    def create_index(self):
        for phrase in self.lang:
            self.vocab.update(phrase.split(' '))

        self.vocab = sorted(self.vocab)

        self.word2idx['<pad>'] = 0
        for index, word in enumerate(self.vocab):
            self.word2idx[word] = index + 1

        for word, index in self.word2idx.items():
            self.idx2word[index] = word


def max_length(tensor):
    return max(len(t) for t in tensor)


def load_dataset(infile):
    # creating cleaned input, output pairs
    # if motif split with whitespace
    if 'motif' in infile:
        pairs = create_dataset_motif(infile)
    # if seq split with list()
    else:
        pairs = create_dataset(infile)


    # index language using the class defined above
    inp_lang = LanguageIndex(first for first, second in pairs)
    targ_lang = LanguageIndex(second for first, second in pairs)

    # Vectorize the input and target languages

    # paratope sequences
    input_tensor = [[inp_lang.word2idx[s] for s in first.split(' ')] for first, second in pairs]

    # epitope sequences
    target_tensor = [[targ_lang.word2idx[s] for s in second.split(' ')] for first, second in pairs]

    # Calculate max_length of input and output tensor
    # Here, we'll set those to the longest sentence in the dataset
    max_length_inp, max_length_tar = max_length(input_tensor), max_length(target_tensor)

    # Padding the input and output tensor to the maximum length
    input_tensor = tf.keras.preprocessing.sequence.pad_sequences(input_tensor,
                                                                 maxlen=max_length_inp,
                                                                 padding='post')

    target_tensor = tf.keras.preprocessing.sequence.pad_sequences(target_tensor,
                                                                  maxlen=max_length_tar,
                                                                  padding='post')

    return input_tensor, target_tensor, inp_lang, targ_lang, max_length_inp, max_length_tar


# run stuff
# make_tabsep_paraepi()
# create_dataset('../../dataset/paraepi.tsv')
infile = '../../dataset/paraepi.tsv'
# infile = '../../dataset/motif_epipara.tsv'
input_tensor, target_tensor, inp_lang, targ_lang, max_length_inp, max_length_targ = load_dataset(infile)
print(input_tensor)
print(inp_lang.idx2word)
print(targ_lang.idx2word)
# sys.exit()
# training and validation sets using an 80-20 split
input_tensor_train, input_tensor_val, target_tensor_train, target_tensor_val = train_test_split(input_tensor, target_tensor, test_size=0.2)

# check length
print(len(input_tensor_train), len(target_tensor_train), len(input_tensor_val), len(target_tensor_val))

### params
BUFFER_SIZE = len(input_tensor_train)
BATCH_SIZE = 32
N_BATCH = BUFFER_SIZE//BATCH_SIZE
embedding_dim = 4
units = 4
vocab_inp_size = len(inp_lang.word2idx)
vocab_tar_size = len(targ_lang.word2idx)

dataset = tf.data.Dataset.from_tensor_slices((input_tensor_train, target_tensor_train)).shuffle(BUFFER_SIZE)
dataset = dataset.batch(BATCH_SIZE, drop_remainder=True)
print(dataset)

def gru(units):
  # If GPU avalaible use CuDNNGRU (cuda stuff)
  # else use vanilla GRU.
  if tf.test.is_gpu_available():
    return tf.keras.layers.CuDNNGRU(units,
                                    return_sequences=True,
                                    return_state=True,
                                    recurrent_initializer='glorot_uniform')
  else:
    return tf.keras.layers.GRU(units,
                               return_sequences=True,
                               return_state=True,
                               recurrent_activation='sigmoid',
                               recurrent_initializer='glorot_uniform')


class Encoder(tf.keras.Model):
    def __init__(self, vocab_size, embedding_dim, enc_units, batch_sz):
        super(Encoder, self).__init__()
        self.batch_sz = batch_sz
        self.enc_units = enc_units
        self.embedding = tf.keras.layers.Embedding(vocab_size, embedding_dim)
        self.gru = gru(self.enc_units)

    def call(self, x, hidden):
        x = self.embedding(x)
        output, state = self.gru(x, initial_state=hidden)
        return output, state

    def initialize_hidden_state(self):
        return tf.zeros((self.batch_sz, self.enc_units))


class Decoder(tf.keras.Model):
    def __init__(self, vocab_size, embedding_dim, dec_units, batch_sz):
        super(Decoder, self).__init__()
        self.batch_sz = batch_sz
        self.dec_units = dec_units
        self.embedding = tf.keras.layers.Embedding(vocab_size, embedding_dim)
        self.gru = gru(self.dec_units)
        self.fc = tf.keras.layers.Dense(vocab_size)

        # used for attention
        self.W1 = tf.keras.layers.Dense(self.dec_units)
        self.W2 = tf.keras.layers.Dense(self.dec_units)
        self.V = tf.keras.layers.Dense(1)

    def call(self, x, hidden, enc_output):
        # enc_output shape == (batch_size, max_length, hidden_size)

        # hidden shape == (batch_size, hidden size)
        # hidden_with_time_axis shape == (batch_size, 1, hidden size)
        # we are doing this to perform addition to calculate the score
        hidden_with_time_axis = tf.expand_dims(hidden, 1)

        # score shape == (batch_size, max_length, 1)
        # we get 1 at the last axis because we are applying tanh(FC(EO) + FC(H)) to self.V
        score = self.V(tf.nn.tanh(self.W1(enc_output) + self.W2(hidden_with_time_axis)))

        # attention_weights shape == (batch_size, max_length, 1)
        attention_weights = tf.nn.softmax(score, axis=1)

        # context_vector shape after sum == (batch_size, hidden_size)
        context_vector = attention_weights * enc_output
        context_vector = tf.reduce_sum(context_vector, axis=1)

        # x shape after passing through embedding == (batch_size, 1, embedding_dim)
        x = self.embedding(x)

        # x shape after concatenation == (batch_size, 1, embedding_dim + hidden_size)
        x = tf.concat([tf.expand_dims(context_vector, 1), x], axis=-1)

        # passing the concatenated vector to the GRU
        output, state = self.gru(x)

        # output shape == (batch_size * 1, hidden_size)
        output = tf.reshape(output, (-1, output.shape[2]))

        # output shape == (batch_size * 1, vocab)
        x = self.fc(output)

        return x, state, attention_weights

    def initialize_hidden_state(self):
        return tf.zeros((self.batch_sz, self.dec_units))

encoder = Encoder(vocab_inp_size, embedding_dim, units, BATCH_SIZE)
decoder = Decoder(vocab_tar_size, embedding_dim, units, BATCH_SIZE)


optimizer = tf.train.AdamOptimizer()


def loss_function(real, pred):
  mask = 1 - np.equal(real, 0)
  loss_ = tf.nn.sparse_softmax_cross_entropy_with_logits(labels=real, logits=pred) * mask
  return tf.reduce_mean(loss_)

checkpoint_dir = './training_checkpoints'
checkpoint_prefix = os.path.join(checkpoint_dir, "ckpt")
checkpoint = tf.train.Checkpoint(optimizer=optimizer,
                                 encoder=encoder,
                                 decoder=decoder)

# # training, split this into a diff code.
# EPOCHS = 10
#
# for epoch in range(EPOCHS):
#     start = time.time()
#
#     hidden = encoder.initialize_hidden_state()
#     total_loss = 0
#
#     for (batch, (inp, targ)) in enumerate(dataset):
#         loss = 0
#         # print(batch, (inp,targ))
#         # sys.exit()
#
#         with tf.GradientTape() as tape:
#             enc_output, enc_hidden = encoder(inp, hidden)
#
#             dec_hidden = enc_hidden
#
#             dec_input = tf.expand_dims([targ_lang.word2idx['<start>']] * BATCH_SIZE, 1)
#
#             # Teacher forcing - feeding the target as the next input
#             for t in range(1, targ.shape[1]):
#                 # passing enc_output to the decoder
#                 predictions, dec_hidden, _ = decoder(dec_input, dec_hidden, enc_output)
#
#                 loss += loss_function(targ[:, t], predictions)
#
#                 # using teacher forcing
#                 dec_input = tf.expand_dims(targ[:, t], 1)
#
#         batch_loss = (loss / int(targ.shape[1]))
#
#         total_loss += batch_loss
#
#         variables = encoder.variables + decoder.variables
#
#         gradients = tape.gradient(loss, variables)
#
#         optimizer.apply_gradients(zip(gradients, variables))
#
#         if batch % 100 == 0:
#             print('Epoch {} Batch {} Loss {:.4f}'.format(epoch + 1,
#                                                          batch,
#                                                          batch_loss.numpy()))
#     # saving (checkpoint) the model every 2 epochs
#     if (epoch + 1) % 2 == 0:
#         checkpoint.save(file_prefix=checkpoint_prefix)
#
#     print('Epoch {} Loss {:.4f}'.format(epoch + 1,
#                                         total_loss / N_BATCH))
#     print('Time taken for 1 epoch {} sec\n'.format(time.time() - start))


def evaluate(sentence, encoder, decoder, inp_lang, targ_lang, max_length_inp, max_length_targ):
    attention_plot = np.zeros((max_length_targ, max_length_inp))

    sentence = preprocess_sentence(sentence)
    # sentence = preprocess_sentence_motif(sentence)

    inputs = [inp_lang.word2idx[i] for i in sentence.split(' ')]
    inputs = tf.keras.preprocessing.sequence.pad_sequences([inputs], maxlen=max_length_inp, padding='post')
    inputs = tf.convert_to_tensor(inputs)

    result = ''

    hidden = [tf.zeros((1, units))]
    enc_out, enc_hidden = encoder(inputs, hidden)

    dec_hidden = enc_hidden
    dec_input = tf.expand_dims([targ_lang.word2idx['<start>']], 0)

    for t in range(max_length_targ):
        predictions, dec_hidden, attention_weights = decoder(dec_input, dec_hidden, enc_out)

        # storing the attention weights to plot later on
        attention_weights = tf.reshape(attention_weights, (-1,))
        attention_plot[t] = attention_weights.numpy()

        predicted_id = tf.argmax(predictions[0]).numpy()

        result += targ_lang.idx2word[predicted_id] + ' '

        if targ_lang.idx2word[predicted_id] == '<end>':
            return result, sentence, attention_plot

        # the predicted ID is fed back into the model
        dec_input = tf.expand_dims([predicted_id], 0)

    return result, sentence, attention_plot


# function for plotting the attention weights
def plot_attention(attention, sentence, predicted_sentence):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1)
    ax.matshow(attention, cmap='viridis')

    fontdict = {'fontsize': 14}
    print(sentence)
    # ax.set_xticklabels([''] + sentence, fontdict=fontdict, rotation=0)
    ax.set_yticklabels([''] + predicted_sentence, fontdict=fontdict)
    plt.xticks(range(0,len(sentence)), list(sentence))
    plt.yticks(range(0,len(predicted_sentence)), list(predicted_sentence))

    plt.show()


def translate(sentence, encoder, decoder, inp_lang, targ_lang, max_length_inp, max_length_targ):
    result, sentence, attention_plot = evaluate(sentence, encoder, decoder, inp_lang, targ_lang, max_length_inp,
                                                max_length_targ)

    print('Input: {}'.format(sentence))
    print('Predicted translation: {}'.format(result))

    attention_plot = attention_plot[:len(result.split(' ')), :len(sentence.split(' '))]
    plot_attention(attention_plot, sentence.split(' '), result.split(' '))


# restoring the latest checkpoint in checkpoint_dir
checkpoint.restore(tf.train.latest_checkpoint(checkpoint_dir))

translate('RAHMADAKAR', encoder, decoder, inp_lang, targ_lang, max_length_inp, max_length_targ)
