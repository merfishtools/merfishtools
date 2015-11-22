import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns


def plot_cdf(values, probs):
    plt.step(values, np.cumsum(probs), "k-")


plt.figure()
sns.set_style("ticks")
pmf = pd.read_table(snakemake.input.pmf, index_col=0)
pmf = pmf.loc[snakemake.wildcards.gene]
pmf.sort("log2fc", inplace=True)

plot_cdf(pmf["log2fc"], np.exp(pmf["prob"]))
ylim = plt.ylim()
plt.xlabel("log2 fold change")
plt.ylabel("CDF")

est = pd.read_table(snakemake.input.est, index_col=0)
est = est.loc[snakemake.wildcards.gene]
ev = est["log2fc_ev"]
sd = est["log2fc_sd"]
plt.vlines([ev], *ylim, colors="red", linestyles="-", label="conditional expectation")
plt.vlines([ev - sd, ev + sd], *ylim, colors="black", linestyles="--", label="standard deviation")
plt.legend(loc="lower right")
plt.text(plt.xlim()[0] + 0.001, ylim[1], "PEP = {:.2}".format(est["diff_pep"]), va="top")

sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")
