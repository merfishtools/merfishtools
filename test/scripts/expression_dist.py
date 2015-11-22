import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns


plt.figure()
sns.set_style("ticks")
exprs = pd.read_table(snakemake.input[0], index_col=0, header=[0, 1])

for cell in exprs:
    expr = exprs[cell]
    expr = np.log10(expr + 1)
    sns.distplot(expr, hist=False, color="black", kde_kws={"alpha": 0.1, "lw": 1, "clip": [0, expr.max()]})

plt.xlim([0, plt.xlim()[1]])
plt.xlabel("log10 expression")
plt.ylabel("density")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
