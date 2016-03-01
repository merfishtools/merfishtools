import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

ests = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input],
                 keys=["{} vs {}".format(*c) for c in snakemake.params.comparisons])

significant = (ests["diff_pep"] <= 0.05).unstack(0).sum(axis="columns")
significant = significant[~(significant.index.str.startswith("notarget") | significant.index.str.startswith("blank"))]
recurrent = significant >= 8
significant[recurrent].to_csv(snakemake.output.foreground, sep="\t")
significant[~recurrent].to_csv(snakemake.output.background, sep="\t")

matrix = ests["log2fc_ev"].unstack(0)
matrix = matrix[~(matrix.index.str.startswith("notarget") | matrix.index.str.startswith("blank"))]

matrix.abs().sum(axis="columns").sort_values(ascending=False).to_csv(snakemake.output.ranked, sep="\t")

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])
cg = sns.clustermap(matrix)
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
plt.savefig(snakemake.output.clust, bbox_inches="tight")
