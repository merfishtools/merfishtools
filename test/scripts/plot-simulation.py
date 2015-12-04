import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])

sim_counts = pd.read_table(snakemake.input.sim_counts, index_col=[0, 1])
known_counts = pd.read_table(snakemake.input.known_counts, index_col=[0, 1])

counts = sim_counts.join(known_counts)
counts["total"] = counts["exact"] + counts["corrected"]
counts["exact_ratio"] = counts["exact"] / counts["total"]

# cmap = mpl.cm.get_cmap("Greys")
# norm = mpl.colors.Normalize(vmin=0, vmax=1)
# plt.subplot(211, aspect="equal")
bins = np.linspace(0, 1, 5)
counts["idx"] = np.digitize(counts["exact_ratio"], bins)
for i, counts in counts.groupby("idx"):
    vmin = bins[i - 1]
    print(counts)
    sns.distplot((counts["total"] + 1) / (counts["known_count"] + 1), hist=False, kde=True, label="≥{}".format(vmin))

plt.xlabel("calling rate")
plt.ylabel("density")
sns.despine()
plt.legend(loc="upper right", title="% exact")
# ax = plt.subplot(212)
# cb = mpl.colorbar.ColorbarBase(ax=ax, cmap=cmap, norm=norm, orientation="vertical")
# cb.set_label("%% exact")

plt.savefig(snakemake.output[0], bbox_inches="tight")
