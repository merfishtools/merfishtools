import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

epsilon = 0.001

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)

plt.figure(figsize=snakemake.config["plots"]["figsize"])
plt.subplot(111, aspect="equal")

errors = []
for i, (mean, posterior_counts, raw_counts, known_counts) in enumerate(zip(
    snakemake.params.means,
    snakemake.input.posterior_counts,
    snakemake.input.raw_counts,
    snakemake.input.known_counts)):

    posterior_counts = pd.read_table(posterior_counts, index_col=[0, 1])["expr_ev"]
    raw_counts = pd.read_table(raw_counts, index_col=[0, 1])
    raw_counts = raw_counts["exact"] + raw_counts["corrected"]
    known_counts = pd.read_table(known_counts, index_col=[0, 1])
    codebook = "mhd4" if snakemake.wildcards.dist == "4" else "mhd2"
    known_counts = known_counts[known_counts[codebook]]

    raw_counts = raw_counts.reindex(known_counts.index, fill_value=0)
    posterior_counts = posterior_counts.reindex(known_counts.index, fill_value=0)

    plt.plot(known_counts["count"] + epsilon, posterior_counts + epsilon, "r.", label="conditional expectation" if i == 0 else "", zorder=1, alpha=0.2, rasterized=True)
    plt.plot(known_counts["count"] + epsilon, raw_counts + epsilon, "k.", label="raw counts" if i == 0 else "", zorder=0, alpha=0.2, rasterized=True)

    errors.append(pd.DataFrame({"error": raw_counts - known_counts["count"], "mean": mean, "type": "raw"}))
    errors.append(pd.DataFrame({"error": posterior_counts - known_counts["count"], "mean": mean, "type": "posterior"}))

plt.plot([epsilon, plt.xlim()[1]], [epsilon, plt.ylim()[1]], "--k")
maxv = max(plt.xlim()[1], plt.ylim()[1])
plt.xlim((0, maxv))
plt.ylim((0,maxv))
plt.ylabel("predicted")
plt.xlabel("truth")
plt.legend(loc="upper left")
sns.despine()
plt.savefig(snakemake.output.scatter, bbox_inches="tight")

errors = pd.concat(errors)
s = (errors["type"] == "posterior") & (errors["mean"] == 5)
print(errors[s].describe())

plt.figure(figsize=snakemake.config["plots"]["figsize"])
colors = sns.xkcd_palette(["light grey", "grey"])
sns.violinplot(x="mean", y="error", hue="type", data=errors, bw=1, split=True, inner="quartile", linewidth=1, palette=colors)
plt.plot(plt.xlim(), [0, 0], "-k", linewidth=1, zorder=-5)

plt.xlabel("mean expression")
plt.ylabel("predicted - truth")
plt.legend(loc="lower left")
sns.despine()

plt.savefig(snakemake.output.violin, bbox_inches="tight")
