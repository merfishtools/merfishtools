import pandas as pd

# load data
exprs = pd.read_table(snakemake.input[0], index_col=[0, 1, 2])["expr_ev"].unstack([0, 1]).fillna(0)
exprs.to_csv(snakemake.output[0], sep="\t")
