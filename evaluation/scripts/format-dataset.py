import numpy as np
import pandas as pd

dataset = pd.read_table(snakemake.input[0], index_col=0)

dataset = dataset.loc[int(snakemake.wildcards.experiment)]

group = snakemake.config["groups"][snakemake.wildcards.group]
dataset = dataset[dataset["Cell_ID"].astype(np.str).str.match(group)]

dataset = dataset.loc[:, ["Cell_ID", "Gene_Name", "Exact_Match", "Cell_Position_X", "Cell_Position_Y", "RNA_Position_X", "RNA_Position_Y"]]
dataset.columns = ["cell", "feat", "dist", "cell_x", "cell_y", "x", "y"]
dataset["dist"] = 1 - dataset["dist"]
dataset.to_csv(snakemake.output[0], sep="\t", index=False)
