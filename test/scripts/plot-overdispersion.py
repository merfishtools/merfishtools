import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])

exprs = pd.read_table(snakemake.input[0], index_col=0, header=[0, 1])

for expmnt, expr in exprs.groupby(level="expmnt", axis=1):
    mean = expr.mean(axis="columns")
    fano = expr.var(axis="columns") / mean
    plt.loglog(mean, fano, "k.", alpha=0.5)
x = plt.xlim()
plt.loglog(x, [1, 1], "-r", alpha=0.8)

plt.xlabel("mean expression")
plt.ylabel("fano factor")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
