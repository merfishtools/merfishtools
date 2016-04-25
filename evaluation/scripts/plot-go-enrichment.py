import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
figsize = snakemake.config["plots"]["figsize"]
figsize[0] *= 3
plt.figure(figsize=figsize)

enrichment = pd.read_table(snakemake.input[0])
enrichment = enrichment[enrichment["adjPvalue"] <= 0.05]

pos = enrichment.index + 0.5

pal = sns.color_palette("Dark2")

plt.barh(pos, enrichment["ExpCount"], align="center", label="expected count", linewidth=0, color=pal[0])
plt.barh(pos, enrichment["Count"], left=enrichment["ExpCount"], align="center", label="observed count", linewidth=0, color=pal[1])
plt.barh(pos, enrichment["Size"], left=enrichment["Count"], align="center", label="term size", linewidth=0, color=pal[2])
plt.yticks(pos, enrichment["Term"])
plt.legend(loc="best")

sns.despine(left=True)
plt.gca().yaxis.set_ticks_position("none")
plt.savefig(snakemake.output[0], bbox_inches="tight")
