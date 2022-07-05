import random
import argparse
import tensorflow
import kerasAC
import pysam
import pandas as pd
import numpy as np
from scipy.special import softmax
from keras.models import load_model
from keras.utils.generic_utils import get_custom_objects
from kerasAC.metrics import *
from kerasAC.custom_losses import *
from tensorflow.keras.models import load_model, model_from_json

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--model_path")
    parser.add_argument("--table")
    parser.add_argument("--fasta")
    parser.add_argument("--out_prefix")
    return parser.parse_args()

def one_hot_encode(seqs):
    """
    Converts a list of DNA ("ACGT") sequences to one-hot encodings, where the
    position of 1s is ordered alphabetically by "ACGT". `seqs` must be a list
    of N strings, where every string is the same length L. Returns an N x L x 4
    NumPy array of one-hot encodings, in the same order as the input sequences.
    All bases will be converted to upper-case prior to performing the encoding.
    Any bases that are not "ACGT" will be given an encoding of all 0s.
    """
    seq_len = len(seqs[0])
    assert np.all(np.array([len(s) for s in seqs]) == seq_len)

    # Join all sequences together into one long string, all uppercase
    seq_concat = "".join(seqs).upper()

    one_hot_map = np.identity(5)[:, :-1]

    # Convert string into array of ASCII character codes;
    base_vals = np.frombuffer(bytearray(seq_concat, "utf8"), dtype=np.int8)

    # Anything that's not an A, C, G, or T gets assigned a higher code
    base_vals[~np.isin(base_vals, np.array([65, 67, 71, 84]))] = 85

    # Convert the codes into indices in [0, 4], in ascending order by code
    _, base_inds = np.unique(base_vals, return_inverse=True)

    # Get the one-hot encoding for those indices, and reshape back to separate
    return one_hot_map[base_inds].reshape((len(seqs), seq_len, 4))

def dinuc_shuffle(seq):
    #get list of dinucleotides
    nucs=[]
    for i in range(0,len(seq),2):
        nucs.append(seq[i:i+2])
    #generate a random permutation
    #set the seed so this shuffling is reproducible for a given sequence
    random.seed(1234)
    random.shuffle(nucs)
    return ''.join(nucs)

def get_preds(model,seq_onehot):
    #seq_onehot = np.expand_dims(seq_onehot, axis=0)
    preds=model.predict(seq_onehot, batch_size=64)
    prof=preds[0] #logits
    probs=np.squeeze(softmax(prof,axis=1)) #probabilities
    count=np.squeeze(preds[1])   #count (1 val)
    count_track=probs*np.exp(count)[:,np.newaxis]  #count track
    return count_track
    #return np.exp(count)

def main():
    args = parse_args()
    model = model_from_json(open(args.model_path + '.arch', 'r').read())
    model.load_weights(args.model_path + '.weights')
    print("loaded_model")

    df = pd.read_csv(args.table, sep='\t')

    chroms = ['chr' + str(x) for x in range(1, 20)] + ['chrX']
    for chrom in chroms:
        print(chrom)
        chrom_df = df.loc[df['chr'] == chrom].reset_index()
        print("Length:", len(chrom_df))

        orig_counts = []

        ref = pysam.FastaFile(args.fasta)

        for i in range(0, len(chrom_df), 10000):
            sub_df = chrom_df.iloc[i:i+10000,:]
    
            orig_onehot = []

            for index,row in sub_df.iterrows():
                center = ((int(row['end']) - int(row['start'])) // 2) + int(row['start'])
                start_index = 1057 - (center - row['start'])
                end_index = 1057 + (row['end'] - center)
                seq = ref.fetch(row['chr'], center - 1057, center + 1057)
                onehot_seq = one_hot_encode(seq)

                orig_onehot.append(onehot_seq)

                if index % 10000 == 0:
                    print(index)
            
            orig_onehot = np.squeeze(np.array(orig_onehot))
            print(orig_onehot.shape)

            orig_counts_loc = get_preds(model, orig_onehot)
            if len(orig_counts) == 0:
                orig_counts = np.array(orig_counts_loc)
            else:
                orig_counts = np.concatenate((orig_counts, orig_counts_loc))

        np.save(args.out_prefix + '/' + chrom + '.orig_counts.npy', orig_counts)

if __name__=="__main__":
    main()

