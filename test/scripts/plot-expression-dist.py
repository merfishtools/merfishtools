import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])
exprs = pd.read_table(snakemake.input[0], index_col=0, header=[0, 1])

for cell in exprs:
    expr = exprs[cell]
    expr = np.log10(expr + 1)
    sns.distplot(expr, hist=False, color="black", kde_kws={"alpha": 0.1, "lw": 1, "clip": [0, expr.max()]})

means = exprs.mean(axis="rows")
sns.distplot(np.log10(means + 1), hist=False, color="red", kde_kws={"clip": [0, means.max()]}, label="cell means")

plt.xlim([0, plt.xlim()[1]])
plt.xlabel("log10 expression")
plt.ylabel("density")
plt.legend(loc="upper right")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
