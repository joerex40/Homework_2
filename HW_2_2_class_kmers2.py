from Bio import SeqIO    # parsing fasta with sequence
import time    # timing
import pandas as pd    # dataframe


class KMer:
    def __init__(self, kmer_name):
        self.sequence = kmer_name
        self.coordinates = []
        self.counter = 1

    def increase(self):
        self.counter += 1

    def add_locus(self, locus):
        self.coordinates.append(locus)

    def data_frame(self):
        df = pd.DataFrame({'start': self.coordinates})
        df['end'] = df.start + len(self.sequence)
        return df


k_mer_size = int(input('Type k-mer length: '))   # determine k-mer length
if k_mer_size == 0:
    print('K-mer length will changed to 1 by default!', '\n')
    k_mer_size = 1


def get_all_k_mer(text, klen):
    """ This function returns a dict of K-Mers objects (with typed length) """
    kmer_dict = {}
    str_range = len(text) - klen + 1
    for i in range(str_range):
        current_kmer = text[i:i + klen]
        if current_kmer in kmer_dict:
            kmer_dict[current_kmer].increase()
            kmer_dict[current_kmer].add_locus(i)
        else:
            kmer_dict[current_kmer] = KMer(current_kmer)
            kmer_dict[current_kmer].add_locus(i)
    return kmer_dict

seq = SeqIO.parse('seq_y_pestis.fasta', "fasta")  # read sequence
string = str(next(seq).seq.upper())   # de-masking string

start = time.time()

k_mer_dict = get_all_k_mer(string, k_mer_size)

end = time.time()
print('Consumed time (sec):', end - start, '\n')

k_mer_champion = max(k_mer_dict.keys(), key=(lambda key: k_mer_dict[key].counter))  # find max k-mer frequency in the dict

df = k_mer_dict[k_mer_champion].data_frame()

print('The most frequent k-mer is:\n{}\n\nCoordinates:\n{}'.format(k_mer_champion, df))
