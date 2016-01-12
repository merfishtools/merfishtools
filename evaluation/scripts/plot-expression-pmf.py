import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns


def plot_pmf(values, probs, stem=False):
    if stem:
        _, _, baseline = plt.stem(values, probs, markerfmt="ko", basefmt="", linefmt="k-")
        plt.setp(baseline, 'linewidth', 0)
    else:
        plt.plot(values, probs, "ko", label="")


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])

pmf = pd.read_table(snakemake.input.expr, index_col=1)
pmf = pmf.loc[snakemake.wildcards.gene]

raw_counts = pd.read_table(snakemake.input.raw_counts, index_col=1)
raw_counts = raw_counts.loc[snakemake.wildcards.gene]

count_exact = raw_counts["exact"]
count_corrected = raw_counts["corrected"]
count_total = count_exact + count_corrected

ylim = (0, np.exp(pmf["prob"]).max())

plt.vlines([count_total], *ylim, colors="red", linestyles="--", label="total count", zorder=5)
plt.vlines([count_exact], *ylim, colors="black", linestyles="--", label="exact count", zorder=5)
plot_pmf(pmf["expr"], np.exp(pmf["prob"]))
plt.ylim(ylim)
plt.xlabel("expression")
plt.ylabel("PMF")
lf = plt.legend(loc="best", frameon=True).get_frame()
lf.set_facecolor("white")
lf.set_edgecolor("white")
plt.xlim([0, plt.xlim()[1]])

sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")
