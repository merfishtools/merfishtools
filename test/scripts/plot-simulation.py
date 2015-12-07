import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])

posterior_counts = pd.read_table(snakemake.input.posterior_counts, index_col=[0, 1])["expr_ev"]
raw_counts = pd.read_table(snakemake.input.raw_counts, index_col=[0, 1])
raw_counts = raw_counts["exact"] + raw_counts["corrected"]
known_counts = pd.read_table(snakemake.input.known_counts, index_col=[0, 1], squeeze=True)

raw_counts = raw_counts.reindex(known_counts.index, fill_value=0)
posterior_counts = posterior_counts.reindex(known_counts.index, fill_value=0)

raw_se = (raw_counts - known_counts) ** 2
posterior_se = (posterior_counts - known_counts) ** 2

print(raw_se.mean(), posterior_se.mean())

sns.distplot(raw_se, kde=True, hist=False, color="black", label="raw counts")
sns.distplot(posterior_se, kde=True, hist=False, color="red", label="conditional expectation")

sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
