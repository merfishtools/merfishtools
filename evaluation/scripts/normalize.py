import pandas as pd
import numpy as np


exprs = [pd.read_table(f, index_col=0) for f in snakemake.input]

# calculate upper quartiles
quartiles = [e.unstack().quantile(0.75) for e in exprs]
# calculate scale factors for upper quartile normalization
scales = pd.Series([quartiles[0] / q for q in quartiles], index=snakemake.params.experiments, name="scale_factors")
scales.to_csv(snakemake.output[0], sep="\t")
