import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])

raw_rmse = []
posterior_rmse = []
for posterior_counts, raw_counts, known_counts in zip(
    snakemake.input.posterior_counts,
    snakemake.input.raw_counts,
    snakemake.input.known_counts):

    posterior_counts = pd.read_table(posterior_counts, index_col=[0, 1])["expr_ev"]
    print(raw_counts)
    raw_counts = pd.read_table(raw_counts, index_col=[0, 1])
    raw_counts = raw_counts["exact"] + raw_counts["corrected"]
    known_counts = pd.read_table(known_counts, index_col=[0, 1], squeeze=True)

    raw_counts = raw_counts.reindex(known_counts.index, fill_value=0)
    posterior_counts = posterior_counts.reindex(known_counts.index, fill_value=0)

    raw_se = (raw_counts - known_counts) ** 2
    posterior_se = (posterior_counts - known_counts) ** 2

    raw_rmse.append(np.sqrt(raw_se.mean()))
    posterior_rmse.append(np.sqrt(posterior_se.mean()))

plt.plot(snakemake.params.means, raw_rmse, "-ko", label="raw counts")
plt.plot(snakemake.params.means, posterior_rmse, "-ro", label="conditional expectation")

plt.xlabel("mean expression")
plt.ylabel("RMSE")
plt.legend(loc="upper left")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
