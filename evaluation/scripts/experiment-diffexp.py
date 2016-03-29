import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

bf_categories = [-np.inf, 0, 2, 6, 10, np.inf]
bf_labels = [0, 1, 2, 3, 4]
#bf_labels = ["no evidence", "weak", "positive", "strong", "very strong"]

ests = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input],
                 keys=["{} vs {}".format(*c) for c in snakemake.params.comparisons])

bf = ests["diff_2lnbf"]
bf = bf.unstack(0)
bf = bf[(bf >= 2).any(axis="columns")]
bf = bf.fillna(-np.inf)
bf = bf.stack()
bf = pd.cut(bf, bf_categories, labels=bf_labels, right=False, include_lowest=True)

matrix = bf.unstack(1)
bf = bf[~(bf.index.str.startswith("notarget") | bf.index.str.startswith("blank"))]

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=np.asarray(snakemake.config["plots"]["figsize"]) * 3)
cg = sns.heatmap(bf)
#plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
plt.savefig(snakemake.output[0], bbox_inches="tight")
