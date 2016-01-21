import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)

small = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input.small], axis=1)
large = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input.large], axis=1)

small_counts = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input.small_counts], axis=1)
large_counts = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input.large_counts], axis=1)

neighbors = pd.read_table(snakemake.input.neighbors, index_col=0, squeeze=True)

common_genes = small.index.intersection(large.index)
means = lambda estimates: estimates.loc[common_genes].mean(axis="columns")

small = means(small)
large = means(large)
small_counts = means(small_counts)
large_counts = means(large_counts)

neighbors = neighbors[common_genes]

min_value = 0.1#min(small.min(), large.min())
max_value = 100#max(small.max(), large.max())

fig = plt.figure(figsize=snakemake.config["plots"]["figsize"])
ax = fig.add_subplot(111, aspect='equal')

ax.scatter(small, large, c=neighbors, cmap=mpl.cm.get_cmap("Reds"), label="conditional expectation", edgecolors="face", norm=mpl.colors.Normalize(vmin=5, vmax=50))
ax.plot([min_value, max_value], [min_value, max_value], "k--")
ax.scatter(small_counts, large_counts, c=neighbors, cmap=mpl.cm.get_cmap("Greys"), label="raw counts", edgecolors="face", norm=mpl.colors.Normalize(vmin=5, vmax=50))
ax.set_xscale("log")
ax.set_yscale("log")

#sns.jointplot(small, large, kind="reg")
plt.xlabel("MHD4 encoding")
plt.ylabel("MHD2 encoding")
plt.xlim((min_value, max_value))
plt.ylim((min_value, max_value))
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
