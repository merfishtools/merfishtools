import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
idx = pd.IndexSlice

import merfishtools


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])

cdf = merfishtools.read_cdf(snakemake.input.fc).loc[snakemake.wildcards.gene]
est = merfishtools.read_diffexp_estimates(snakemake.input.fc_est).loc[snakemake.wildcards.gene]

merfishtools.plot_cdf(cdf, expected_value=est["log2fc_ev"], credible_interval=est[["log2fc_ci_lower", "log2fc_ci_upper"]], legend=snakemake.wildcards.legend == "legend")

plt.xlabel("log2 fold change")
plt.ylabel("CDF")
if snakemake.wildcards.legend == "legend":
    plt.legend(loc="best")

sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")
