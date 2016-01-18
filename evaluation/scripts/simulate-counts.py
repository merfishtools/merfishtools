from collections import Counter, defaultdict
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
    return (word ^ err10) ^ err01, err10.count(True) + err01.count(True)


def hamming1_env(word):
    for i in range(len(word)):
        w = word.copy()
        w[i] ^= 1
        yield w


codebook_mhd4 = pd.read_table(snakemake.input.mhd4, index_col=0, dtype=np.dtype(str), squeeze=True).apply(bitarray)
codebook_mhd2 = pd.read_table(snakemake.input.mhd2, index_col=0, dtype=np.dtype(str), squeeze=True).apply(bitarray)
genes = set(codebook_mhd2.index) | set(codebook_mhd4.index)


for gene, a in codebook_mhd4.items():
    neighbors = sum((a ^ b).count(True) == 4 for _, b in codebook_mhd4.items())
    print(gene, neighbors)


with open(snakemake.output.known_counts, "w") as known_out:
    known_out = csv.writer(known_out, delimiter="\t")
    known_out.writerow(["cell", "feat", "mhd2", "mhd4", "count"])
    known_counts = defaultdict(dict)
    for cell in range(snakemake.params.cell_count):
        random_counts = np.random.poisson(int(snakemake.wildcards.mean), len(genes))
        for gene, count in zip(genes, random_counts):
            known_counts[cell][gene] = count
            known_out.writerow([cell, gene, gene in codebook_mhd2.index, gene in codebook_mhd4.index, count])


def simulate(codebook, counts_path, has_corrected=True):
    lookup_exact = {word.tobytes(): gene for gene, word in codebook.items()}
    if has_corrected:
        lookup_corrected = {w.tobytes(): gene for gene, word in codebook.items() for w in hamming1_env(word)}
    else:
        lookup_corrected = {}

    with open(counts_path, "w") as sim_out:
        sim_out = csv.writer(sim_out, delimiter="\t")
        sim_out.writerow(["cell", "feat", "dist", "cell_x", "cell_y", "x", "y"])

        errors = []
        for cell in range(snakemake.params.cell_count):
            readouts = []

            for gene, word in codebook.items():
                count = known_counts[cell][gene]
                if True: #not gene.startswith("notarget") and not gene.startswith("blank"):
                    for _ in range(count):
                        readout, errs = sim_errors(word)
                        errors.append(errs)
                        readouts.append(readout)

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


simulate(codebook_mhd4, snakemake.output.sim_counts_mhd4, has_corrected=True)
simulate(codebook_mhd2, snakemake.output.sim_counts_mhd2, has_corrected=False)
