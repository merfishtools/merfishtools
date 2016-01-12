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

common_genes = small.index.intersection(large.index)

small = small.loc[common_genes].mean(axis="columns")
large = large.loc[common_genes].mean(axis="columns")
small_counts = small_counts.loc[common_genes].mean(axis="columns")
large_counts = large_counts.loc[common_genes].mean(axis="columns")

min_value = 0.1#min(small.min(), large.min())
max_value = 100#max(small.max(), large.max())

fig = plt.figure(figsize=snakemake.config["plots"]["figsize"])
ax = fig.add_subplot(111, aspect='equal')

ax.loglog(small, large, "r.", label="conditional expectation")
ax.loglog([min_value, max_value], [min_value, max_value], "k--")
ax.loglog(small_counts, large_counts, "k.", label="raw counts")

#sns.jointplot(small, large, kind="reg")
plt.xlabel("MHD4 encoding")
plt.ylabel("MHD2 encoding")
plt.xlim((min_value, max_value))
plt.ylim((min_value, max_value))
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
