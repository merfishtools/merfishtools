import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)

plt.figure(figsize=snakemake.config["plots"]["figsize"])


pmfs = []
for i, (mean, pmf, known_counts) in enumerate(zip(
    snakemake.params.means,
    snakemake.input.pmfs,
    snakemake.input.known_counts)):
    print(mean)

    pmf = pd.read_table(pmf, index_col=[0, 1])
    known_counts = pd.read_table(known_counts, index_col=[0, 1])
    codebook = "mhd4" if snakemake.wildcards.dist == "4" else "mhd2"
    known_counts = known_counts[known_counts[codebook]]["count"]
    common = pmf.index.intersection(known_counts.index)

    known_counts = known_counts.loc[common]
    pmf = pmf.loc[common]

    pmf["expr"] -= known_counts
    pmf.set_index("expr", append=True, inplace=True)
    pmf = pmf.unstack(level=[0, 1])
    pmf = np.exp(pmf.sample(10, axis=1))
    pmf.columns = np.arange(pmf.shape[1])
    pmfs.append(pmf)

print("plotting")
pmfs = pd.concat(pmfs, axis=1)
pmfs = pmfs.reindex(np.arange(-10, 11))
pmfs.fillna(0, inplace=True)
pmfs.sort_index(ascending=False, inplace=True)
sns.heatmap(pmfs, cbar=False, cmap="Greys", xticklabels=False, yticklabels=10)

plt.ylabel("predicted - truth")
plt.xlabel("PMF")

plt.savefig(snakemake.output[0], bbox_inches="tight")
