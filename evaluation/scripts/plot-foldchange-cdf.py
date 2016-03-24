import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])

pmf = pd.read_table(snakemake.input.pmf, index_col=0).loc[snakemake.wildcards.gene]
pmf.sort("log2fc", inplace=True)

cdf = np.cumsum(np.exp(pmf["prob"]))
plt.step(pmf["log2fc"], cdf, "k-", label="", clip_on=False, zorder=6)
ylim = plt.ylim()
plt.xlabel("log2 fold change")
plt.ylabel("CDF")

est = pd.read_table(snakemake.input.est, index_col=0).loc[snakemake.wildcards.gene]
ev = est["log2fc_ev"]
ci_lower, ci_upper = est[["log2fc_ci_lower", "log2fc_ci_upper"]]

plt.fill([ci_lower, ci_upper, ci_upper, ci_lower], [0, 0, ylim[1], ylim[1]], "red", lw=0, label="95% credible interval", alpha=0.5)
plt.vlines([ev], *ylim, colors="red", linestyles="-", label="expected value")

plt.fill_between(pmf["log2fc"], cdf, 1.2,  zorder=5, facecolor="white", edgecolor="white", step="pre")

plt.ylim(ylim)
if snakemake.wildcards.legend == "legend":
    plt.legend(loc="best")

sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")
