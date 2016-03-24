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

pmf = pd.read_table(snakemake.input.expr, index_col=1).loc[snakemake.wildcards.gene]
est = pd.read_table(snakemake.input.expr_est, index_col=1).loc[snakemake.wildcards.gene]
raw_counts = pd.read_table(snakemake.input.raw_counts, index_col=1).loc[snakemake.wildcards.gene]

count_exact = raw_counts["exact"]
count_corrected = raw_counts["corrected"]
count_total = count_exact + count_corrected

plot_pmf(pmf["expr"], np.exp(pmf["prob"]))

ylim = plt.ylim()

ci_lower, ci_upper = est[['expr_ci_lower', 'expr_ci_upper']]
plt.fill([ci_lower, ci_upper, ci_upper, ci_lower], [0, 0, ylim[1], ylim[1]], "red", lw=0, label="95% CI", alpha=0.5)
plt.vlines([est['expr_ev']], *ylim, colors="red", label="expected value")

plt.vlines([count_total], *ylim, colors="grey", linestyles="--", label="total count")
plt.vlines([count_exact], *ylim, colors="grey", linestyles=":", label="exact count")

plt.ylim(ylim)
plt.xlabel("expression")
plt.ylabel("PMF")
if snakemake.wildcards.legend == "legend":
    plt.legend(loc="upper left", bbox_to_anchor=(0.5, 1))
plt.xlim([0, plt.xlim()[1]])

sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")
