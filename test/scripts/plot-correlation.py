import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
import random

random.seed(2139)


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
fig = plt.figure(figsize=snakemake.config["plots"]["figsize"])

ax = fig.add_subplot(111, aspect='equal')
exprs = pd.read_table(snakemake.input[0], index_col=0, header=[0, 1])

expr_a = pd.Series()
expr_b = pd.Series()

for expmnt, expr in exprs.groupby(level="expmnt", axis=1):
    if int(expmnt) == 1:
        for a, b in combinations(expr.columns, 2):
            expr_a = expr_a.append(expr[a])
            expr_b = expr_b.append(expr[b])

expr_a = np.log10(1 + expr_a)
expr_b = np.log10(1 + expr_b)
dropout = np.logical_or(np.logical_and(expr_a > 0.2, expr_b == 0), np.logical_and(expr_a == 0, expr_b > 0.2))
high = np.logical_or(expr_a >= 10, expr_b >= 10)
#sns.kdeplot(expr_a[~dropout], expr_b[~dropout], shade=True, cmap="Greys", shade_lowest=False, ax=ax, clip=[0, 1.5])
#sns.kdeplot(expr_a[dropout], expr_b[dropout], shade=True, cmap="Reds", shade_lowest=False, ax=ax, clip=[0, 1.5])
sns.kdeplot(expr_a, expr_b, shade=True, cmap="Greys", shade_lowest=False, ax=ax)

plt.xlabel("log10 expression")
plt.ylabel("log10 expression")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
