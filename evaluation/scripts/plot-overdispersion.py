import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])
ax = plt.subplot(111, aspect="equal")

for f in snakemake.input:
    exprs = pd.read_table(f, index_col=0)
    mean = exprs.mean(axis="columns")
    var = exprs.var(axis="columns")
    plt.loglog(mean, var, "k.", alpha=0.5)
max_val = max(plt.xlim()[1], plt.ylim()[1])
min_val = min(plt.xlim()[0], plt.ylim()[0])
lim = (min_val, max_val)
plt.xlim(lim)
plt.ylim(lim)

x = plt.xlim()
plt.loglog(x, x, "-r", alpha=0.8)

plt.xlabel("mean expression")
plt.ylabel("variance")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
