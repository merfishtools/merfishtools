import pandas as pd
import numpy as np

def hamming_dist(a, b):
    return sum(x != y for x, y in zip(a, b))

codebook = pd.read_table(snakemake.input[0], dtype=np.dtype(str))

full_codebook = codebook
codebook = codebook[~codebook["Gene"].str.startswith("blank") & ~codebook["Gene"].str.startswith("notarget")]

def count_neighbors(codebook):
    return pd.DataFrame({
        "gene": codebook["Gene"],
        "neighbors": [sum(hamming_dist(a, b) <= snakemake.params.neighbor_dist for b in codebook["Codeword"]) for a in codebook["Codeword"]]
    })

neighbors = count_neighbors(codebook)
print(count_neighbors(full_codebook))
#print(neighbors, len(neighbors), np.bincount(neighbors["neighbors"]))

neighbors.to_csv(snakemake.output[0], sep="\t", index=False)
