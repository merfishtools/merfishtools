from collections import Counter
from itertools import chain
import csv

import pandas as pd
import numpy as np
from bitarray import bitarray


p0 = 0.04
p1 = 0.1
dist = int(snakemake.wildcards.dist)


def sim_errors(word):
    p = np.random.uniform(0, 1, len(word))
    err10 = bitarray(list(p <= p1)) & word
    err01 = bitarray(list(p <= p0)) & ~word
    return (word ^ err10) ^ err01, err10.count(True) + err01.count(True)


def hamming1_env(word):
    for i in range(len(word)):
        w = word.copy()
        w[i] ^= 1
        yield w


codebook = pd.read_table(snakemake.input[0], index_col=0, dtype=np.dtype(str), squeeze=True).apply(bitarray)
lookup_exact = {word.tobytes(): gene for gene, word in codebook.items()}
if dist == 4:
    lookup_corrected = {w.tobytes(): gene for gene, word in codebook.items() for w in hamming1_env(word)}
else:
    lookup_corrected = {}

with open(snakemake.output.known_counts, "w") as known_out, open(snakemake.output.sim_counts, "w") as sim_out:
    known_out = csv.writer(known_out, delimiter="\t")
    known_out.writerow(["cell", "feat", "known_count"])

    sim_out = csv.writer(sim_out, delimiter="\t")
    sim_out.writerow(["cell", "feat", "dist", "cell_x", "cell_y", "x", "y"])

    errors = []
    for cell in range(snakemake.params.cell_count):
        readouts = []
        random_counts = np.random.poisson(int(snakemake.wildcards.mean), len(codebook))

        for (gene, word), count in zip(codebook.items(), random_counts):
            if not gene.startswith("notarget") and not gene.startswith("blank"):
                for _ in range(count):
                    readout, errs = sim_errors(word)
                    errors.append(errs)
                    readouts.append(readout)
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
            for _ in range(exact_counts[gene]):
                sim_out.writerow([cell, gene, 0, 0, 0, 0, 0])
            for _ in range(corrected_counts[gene]):
                sim_out.writerow([cell, gene, 1, 0, 0, 0, 0])
    errors = pd.Series(errors)
    print(errors.value_counts())
    print(errors.value_counts() / errors.sum())
