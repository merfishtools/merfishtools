from collections import Counter
from itertools import chain
import csv

import pandas as pd
import numpy as np
from bitarray import bitarray


p0 = 0.04
p1 = 0.1


def sim_errors(word):
    p = np.random.uniform(0, 1, len(word))
    err10 = bitarray(list(p <= p1)) & word
    err01 = bitarray(list(p <= p0)) & ~word
    return (word ^ err10) ^ err01

def hamming1_env(word):
    for i in range(len(word)):
        w = word.copy()
        w[i] ^= 1
        yield w

codebook = pd.read_table(snakemake.input[0], index_col=0, dtype=np.dtype(str), squeeze=True).apply(bitarray)
lookup_exact = {word.tobytes(): gene for gene, word in codebook.items()}
lookup_corrected = {w.tobytes(): gene for gene, word in codebook.items() for w in hamming1_env(word)}

with open(snakemake.output.known_counts, "w") as known_out, open(snakemake.output.sim_counts, "w") as sim_out:
    known_out = csv.writer(known_out, delimiter="\t")
    known_out.writerow(["cell", "gene", "known_count"])

    sim_out = csv.writer(sim_out, delimiter="\t")
    sim_out.writerow(["cell", "gene", "exact", "corrected"])

    for cell in range(snakemake.params.cell_count):
        readouts = []
        random_counts = np.random.poisson(15, len(codebook))

        for (gene, word), count in zip(codebook.items(), random_counts):
            if not gene.startswith("notarget") and not gene.startswith("blank"):
                readouts.extend(sim_errors(word) for _ in range(count))
            else:
                count = 0
            known_out.writerow([cell, gene, count])

        exact_counts = Counter()
        corrected_counts = Counter()
        for readout in readouts:
            try:
                exact_counts[lookup_exact[readout.tobytes()]] += 1
            except KeyError:
                try:
                    corrected_counts[lookup_corrected[readout.tobytes()]] += 1
                except KeyError:
                    # readout is lost
                    pass

        for gene in set(chain(exact_counts, corrected_counts)):
            sim_out.writerow([cell, gene, exact_counts[gene], corrected_counts[gene]])
