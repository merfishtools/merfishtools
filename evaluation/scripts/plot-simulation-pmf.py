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
    for _, _pmf in pmf.groupby(level=[0, 1]):
        _pmf.reset_index(drop=True, inplace=True)
        pmfs.append(pd.Series(data=_pmf["prob"], index=_pmf["expr"]))
    break
        

pmfs = pd.concat(pmfs, axis=1).fillna(0)
sns.heatmap(pmfs)
plt.ylabel("predicted - truth")
plt.xlabel("PMF")

plt.savefig(snakemake.output[0], bbox_inches="tight")
