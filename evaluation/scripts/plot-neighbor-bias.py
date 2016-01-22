import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

exprs = [pd.read_table(f, index_col=0) for f in snakemake.input.counts]
neighbors = pd.read_table(snakemake.input.neighbors, index_col=0, squeeze=True)

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
ax = plt.subplot(111, aspect="equal")

for expr in exprs:
    for cell in expr:
        ax.scatter(neighbors, expr.loc[neighbors.index, cell], c="k", edgecolors="face", rasterized=True)
plt.xlabel("neighbors")
plt.ylabel("counts")
#ax.set_yscale("log")
sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")
