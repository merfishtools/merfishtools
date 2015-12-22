import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])

errors = []
for mean, posterior_counts, raw_counts, known_counts in zip(
    snakemake.params.means,
    snakemake.input.posterior_counts,
    snakemake.input.raw_counts,
    snakemake.input.known_counts):

    posterior_counts = pd.read_table(posterior_counts, index_col=[0, 1])["expr_ev"]
    raw_counts = pd.read_table(raw_counts, index_col=[0, 1])
    raw_counts = raw_counts["exact"] + raw_counts["corrected"]
    known_counts = pd.read_table(known_counts, index_col=[0, 1])
    codebook = "mhd4" if snakemake.wildcards.dist == "4" else "mhd2"
    known_counts = known_counts[known_counts[codebook]]

    raw_counts = raw_counts.reindex(known_counts.index, fill_value=0)
    posterior_counts = posterior_counts.reindex(known_counts.index, fill_value=0)

    errors.append(pd.DataFrame({"error": raw_counts - known_counts["count"], "mean": mean, "type": "raw"}))
    errors.append(pd.DataFrame({"error": posterior_counts - known_counts["count"], "mean": mean, "type": "posterior"}))

errors = pd.concat(errors)

colors = sns.xkcd_palette(["light grey", "grey"])
sns.violinplot(x="mean", y="error", hue="type", data=errors, split=True, inner="quartile", linewidth=1, palette=colors)
plt.plot(plt.xlim(), [0, 0], "-k", linewidth=1, zorder=-5)

plt.xlabel("mean expression")
plt.ylabel("predicted - true")
plt.legend(loc="lower left")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
