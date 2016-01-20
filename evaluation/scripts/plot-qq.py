import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)


quantiles = np.linspace(0, 1, 100)

exprs = [pd.read_table(f, index_col=0) for f in snakemake.input]

plt.figure(figsize=(20, 20))
for (i, (a, a_name)), (j, (b, b_name)) in combinations(enumerate(zip(exprs, snakemake.params.experiments)), 2):
    plt.subplot2grid([len(exprs)] * 2, (i, j))
    q_a = a.unstack().quantile(quantiles)
    q_b = b.unstack().quantile(quantiles)

    vmax = max(q_a.max(), q_b.max())

    plt.loglog(q_a, q_b, "k.")
    plt.loglog([0.01, vmax], [0.01, vmax], "r--")
    plt.xlabel("experiment {}".format(a_name))
    plt.ylabel("experiment {}".format(b_name))
    sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
