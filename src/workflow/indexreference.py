import pickle
import numpy as np
from Bio import SeqIO


def index_reference(ref: str, output: str):
    """
    Generate Reference Dictionary for A/T locations
    :param ref: Name of the Reference File (FASTA)
    :param output: Name of the Index File (pickle)
    :return:
    """

    assert isinstance(ref, str)
    assert isinstance(output, str)

    ref_array = {}
    for record in SeqIO.parse(ref, 'fasta'):
        name = record.id
        seq = str(record.seq).upper()
        seqlen = len(seq)
        tmp = np.zeros(seqlen, dtype=np.bool)
        for i, nt in enumerate(seq):
            if nt == 'A' or nt == 'T':
                tmp[i] = 1
        ref_array[name] = tmp
    with open(output, 'wb') as filep:
        pickle.dump(ref_array, filep, protocol=pickle.HIGHEST_PROTOCOL)


