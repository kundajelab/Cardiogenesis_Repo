import argparse
import pyBigWig
import numpy as np
import h5py
import pickle as pkl
import pysam
import pandas as pd
import pickle

# need full paths!
parser = argparse.ArgumentParser(description="Convert importance scores in numpy format to bigwig. The output can be visualised using WashU Epigenome Browser as a dynseq track. Please read all parameter argument requirements. PROVIDE ABSOLUTE PATHS!")
parser.add_argument("-bed",type=str, required=True, help="Peaks file")
parser.add_argument("-npy",type=str, required=True, help="Numpy file with importance scores")
parser.add_argument("-fasta",type=str, required=True, help="Reference genome fasta")
parser.add_argument("-c", "--chrom_sizes", type=str, required=True, help="Chromosome sizes 2 column tab-separated file")
parser.add_argument("-o", "--outfile", type=str, required=True, help="Output bigwig file")
parser.add_argument("-s", "--outstats", type=str, required=True, help="Output file with stats of low and high quantiles")
parser.add_argument("-t", "--tqdm", type=int,default=0, help="Use tqdm. If yes then you need to have it installed.")
parser.add_argument("--type",type=str)
args = parser.parse_args()

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

print(args)

ref = pysam.FastaFile(args.fasta)

with open(args.chrom_sizes) as f:
    gs = [x.strip().split('\t') for x in f if len(x.strip())!=0]

gs = [(x[0], int(x[1])) for x in gs]

chr_to_idx = {}
for i,x in enumerate(gs):
    chr_to_idx[x[0]] = i

peaks = pd.read_csv(args.bed, sep='\t', header=None)
shaps = np.load(args.npy)

f =  {'seq':{}, args.type:{}}

for index,row in peaks.iterrows():
    chrom = row[0]
    summit = row[1] + row[9]
    start = summit - 250
    end = summit + 250
    seq = ref.fetch(chrom, start, end)
    onehot_seq = np.squeeze(one_hot_encode(seq))
    shap_seq = shaps[index][250:750]
    seq_tup = (chrom, summit)
    f['seq'][seq_tup] = onehot_seq
    f[args.type][seq_tup] = shap_seq

    if index % 10000 == 0:
        print(index)

#with open('/home/soumyak/cad/bpnet_models/baseline_bpnet/test.pkl', 'wb') as outfile:
#    pickle.dump(f, outfile)

SEQLEN = np.array(list(f[args.type].values())).shape[1]
print(SEQLEN)
assert(SEQLEN%2==0)

#with open(args.regions) as r:
#    regions = [x.strip().split('\t') for x in r]
regions = [key for key in f[args.type].keys()]

#regions = [[x[0], int(x[1])-int(SEQLEN/2), int(x[1])+int(SEQLEN/2)] for x in regions]
regions = [[x[0], int(x[1])-250, int(x[1])+250] for x in regions]

# regions may not be sorted, so get their sorted order
order_of_regs = sorted(range(len(regions)), key=lambda x:(chr_to_idx[regions[x][0]], regions[x][1]))

# regions may overlap but as we go in sorted order, we will ignore the values that are repeated
# and only consider the first instance

bw = pyBigWig.open(args.outfile, 'w')
bw.addHeader(gs)

all_entries = []
cur_chr = ""
cur_end = 0

if args.tqdm:
    from tqdm import tqdm
    iterator = tqdm(order_of_regs)
else:
    iterator = order_of_regs

for i in iterator:
    # subset to chromosome (debugging)
    #if regions[i][0]!="chr12":
    #    continue

    if regions[i][0]!=cur_chr:
        cur_chr = regions[i][0]
        cur_end = 0

    # bring current end to at least start of current region
    if cur_end < regions[i][1]:
        cur_end = regions[i][1]

    assert(regions[i][2]>=cur_end)

    tup1 = (cur_chr, int(regions[i][1]+regions[i][2])/2)
    if "shap" in args.type:
        vals = np.array((f[args.type][tup1]), dtype=np.float64)
    else:
        vals = np.array((f[args.type][tup1]*np.exp(f[args.type.replace("prof","sum")][tup1])), dtype=np.float64)
    #print(vals)
    #print(vals.shape)
    #print(len([regions[i][0]]*(regions[i][2]-cur_end)))
    #print(len(list(range(cur_end,regions[i][2]))))
    #print(list(vals))
    bw.addEntries([regions[i][0]]*(regions[i][2]-regions[i][1]),
                   list(range(regions[i][1],regions[i][2])),
                   ends =list(range(regions[i][1]+1, regions[i][2]+1)),
                   values=list(vals))

    all_entries.append(vals)
    cur_end = regions[i][2]+1

bw.close()
#f.close()

all_entries = np.hstack(all_entries)

with open(args.outstats, 'w') as f:
    f.write("Min\t{:.6f}\n".format(np.min(all_entries)))
    f.write(".1%\t{:.6f}\n".format(np.quantile(all_entries, 0.001)))
    f.write("1%\t{:.6f}\n".format(np.quantile(all_entries, 0.01)))
    f.write("50%\t{:.6f}\n".format(np.quantile(all_entries, 0.5)))
    f.write("99%\t{:.6f}\n".format(np.quantile(all_entries, 0.99)))
    f.write("99.9%\t{:.6f}\n".format(np.quantile(all_entries, 0.999)))
    f.write("99.95%\t{:.6f}\n".format(np.quantile(all_entries, 0.9995)))
    f.write("99.99%\t{:.6f}\n".format(np.quantile(all_entries, 0.9999)))
    f.write("Max\t{:.6f}\n".format(np.max(all_entries)))

