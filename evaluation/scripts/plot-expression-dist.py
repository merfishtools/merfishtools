import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])

means = []
for f in snakemake.input:
    exprs = pd.read_table(f, index_col=0)
    for cell in exprs:
        expr = exprs[cell]
        means.append(expr.mean())
        expr = np.log10(expr + 1)
        sns.kdeplot(expr, color="black", clip=[0, expr.max()], alpha=0.1, lw=1, label="", kernel="gau", bw="scott")

means = np.array(means)
sns.kdeplot(np.log10(means + 1), color="red", clip=[0, means.max()], label="cell means")

plt.xlim([0, plt.xlim()[1]])
plt.xlabel("log10 expression")
plt.ylabel("density")
plt.legend(loc="upper right")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
