import pandas as pd
import numpy as np

raw = pd.read_table(snakemake.input[0], index_col=[0, 1])

def exact(d):
    return (1 - d).sum()

def corrected(d):
    return d.sum()

genes = raw.groupby(level=[0, 1])
counts = pd.DataFrame({"exact": genes["dist"].aggregate(exact),
                       "corrected": genes["dist"].aggregate(corrected)})

counts.to_csv(snakemake.output[0], sep="\t")
