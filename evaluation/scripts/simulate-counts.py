from collections import Counter, defaultdict
from itertools import chain
import csv

import pandas as pd
import numpy as np
from bitarray import bitarray


# Fix the seed for a particular mean, so that results are reproducible.
np.random.seed(int(snakemake.wildcards.mean))

p0 = 0.04
p1 = 0.1


def sim_errors(word):
    p = np.random.uniform(0, 1, len(word))
    err10 = bitarray(list(p <= p1)) & word
    err01 = bitarray(list(p <= p0)) & ~word
    readout = (word ^ err10) ^ err01
    errs = err10.count(True) + err01.count(True)
    assert errs == 0 or word != readout
    return readout, errs


def hamming1_env(word):
    for i in range(len(word)):
        w = word.copy()
        w[i] ^= 1
        yield w


codebook_mhd4 = pd.read_table(snakemake.input.mhd4, index_col=0, dtype=np.dtype(str))["codeword"].apply(bitarray)
codebook_mhd2 = pd.read_table(snakemake.input.mhd2, index_col=0, dtype=np.dtype(str))["codeword"].apply(bitarray)
genes = set(codebook_mhd2.index) | set(codebook_mhd4.index)


for gene, a in codebook_mhd4.items():
    neighbors = sum((a ^ b).count(True) == 4 for _, b in codebook_mhd4.items())


with open(snakemake.output.known_counts, "w") as known_out:
    known_out = csv.writer(known_out, delimiter="\t")
    known_out.writerow(["cell", "feat", "mhd2", "mhd4", "count"])
    known_counts = defaultdict(lambda: defaultdict(int))
    for cell in range(snakemake.params.cell_count):
        random_counts = np.random.poisson(int(snakemake.wildcards.mean), len(genes))
        for gene, count in zip(genes, random_counts):
            if True: # not gene.startswith("notarget") and not gene.startswith("blank"):
                known_counts[cell][gene] = count
                known_out.writerow([cell, gene, gene in codebook_mhd2.index, gene in codebook_mhd4.index, count])


def simulate(codebook, counts_path, stats_path, has_corrected=True):
    lookup_exact = {word.tobytes(): gene for gene, word in codebook.items()}
    if has_corrected:
        lookup_corrected = {w.tobytes(): gene for gene, word in codebook.items() for w in hamming1_env(word)}
    else:
        lookup_corrected = {}

    with open(counts_path, "w") as sim_out:
        sim_out = csv.writer(sim_out, delimiter="\t")
        sim_out.writerow(["cell", "feat", "dist", "cell_x", "cell_y", "x", "y"])

        stats = []
        for cell in range(snakemake.params.cell_count):
            readouts = []
            words = []
            genes = []
            errors = []

            for gene, word in codebook.items():
                count = known_counts[cell][gene]
                for _ in range(count):
                    readout, errs = sim_errors(word)
                    errors.append(errs)
                    readouts.append(readout)
                    words.append(word)
                    genes.append(gene)

            exact_counts = Counter()
            corrected_counts = Counter()
            exact_miscalls = Counter()
            corrected_miscalls = Counter()
            for readout, errs, word, orig_gene in zip(readouts, errors, words, genes):
                try:
                    gene = lookup_exact[readout.tobytes()]
                    exact_counts[gene] += 1
                    if errs > 0:
                        exact_miscalls[gene] += 1
                except KeyError:
                    try:
                        corrected_counts[lookup_corrected[readout.tobytes()]] += 1
                        if errs > 1:
                            corrected_miscalls[gene] += 1
                    except KeyError:
                        # readout is lost
                        pass

            for gene in set(chain(exact_counts, corrected_counts)):
                for _ in range(exact_counts[gene]):
                    sim_out.writerow([cell, gene, 0, 0, 0, 0, 0])
                for _ in range(corrected_counts[gene]):
                    sim_out.writerow([cell, gene, 1, 0, 0, 0, 0])

            for gene in exact_counts:
                known = known_counts[cell][gene]
                counts = exact_counts[gene] + corrected_counts[gene]
                _exact_miscalls = exact_miscalls[gene]
                _corrected_miscalls = corrected_miscalls[gene]
                miscalls = _exact_miscalls + _corrected_miscalls
                missed = known - counts + miscalls
                stats.append([cell, gene, known, missed, counts, miscalls])

        stats = pd.DataFrame(stats)
        stats.columns = ["cell", "gene", "truth", "missed", "counts", "miscalls"]
        stats["total"] = stats["truth"] + stats["miscalls"]
        stats["calls"] = stats["counts"] - stats["miscalls"]
        stats["call-rate"] = stats["calls"] / stats["total"]
        stats["miscall-rate"] = stats["calls"] / stats["total"]
        stats["missed-rate"] = stats["missed"] / stats["total"]
        stats.to_csv(stats_path, sep="\t")

        print("rates vs counts")
        print("miscalls", (stats["miscalls"] / stats["counts"]).mean())
        print("calls", (stats["calls"] / stats["counts"]).mean())
        print("missed", (stats["missed"] / stats["counts"]).mean())
        print("rates vs truth")
        print("miscalls", (stats["miscalls"] / stats["truth"]).mean())
        print("calls", (stats["calls"] / stats["truth"]).mean())
        print("missed", (stats["missed"] / stats["truth"]).mean())
        print("rates vs total")
        print("miscalls", (stats["miscalls"] / stats["total"]).mean())
        print("calls", (stats["calls"] / stats["total"]).mean())
        print("missed", (stats["missed"] / stats["total"]).mean())



print("Simulating MHD4")
simulate(codebook_mhd4, snakemake.output.sim_counts_mhd4, snakemake.output.stats_mhd4, has_corrected=True)
print("Simulating MHD2")
simulate(codebook_mhd2, snakemake.output.sim_counts_mhd2, snakemake.output.stats_mhd2, has_corrected=False)
