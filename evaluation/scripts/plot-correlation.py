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

expr_a = pd.Series()
expr_b = pd.Series()
for f in snakemake.input[:1]:
    exprs = pd.read_table(f, index_col=0)
    for a, b in combinations(exprs.columns, 2):
        expr_a = expr_a.append(exprs[a])
        expr_b = expr_b.append(exprs[b])

expr_a = np.log10(1 + expr_a)
expr_b = np.log10(1 + expr_b)
dropout = np.logical_or(np.logical_and(expr_a > 0.2, expr_b <= 0.2), np.logical_and(expr_a <= 0.2, expr_b > 0.2))
high = np.logical_or(expr_a >= 10, expr_b >= 10)
sns.kdeplot(expr_a[~dropout], expr_b[~dropout], shade=True, cmap="Greys", shade_lowest=False, ax=ax, clip=[0, 2])
sns.kdeplot(expr_a[dropout], expr_b[dropout], shade=True, cmap="Reds", shade_lowest=False, ax=ax, clip=[0, 2])
#sns.kdeplot(expr_a, expr_b, shade=True, cmap="Greys", shade_lowest=False, ax=ax)

plt.xlabel("log10 expression")
plt.ylabel("log10 expression")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
