import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns


def hamming_dist(a, b):
    return sum(x != y for x, y in zip(a, b))


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])

codebook = pd.read_table(snakemake.input[0], index_col=0, dtype=np.dtype(str))
neighbors = [sum(hamming_dist(a, b) == snakemake.params.neighbor_dist for b in codebook["Codeword"]) for a in codebook["Codeword"]]
print(neighbors, len(neighbors))

sns.distplot(neighbors, hist=True, kde=False)

plt.savefig(snakemake.output[0], bbox_inches="tight")
