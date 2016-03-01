import pandas as pd
from scipy.spatial import ConvexHull

data = pd.read_table(snakemake.input[0], index_col=0)

def get_area(rnas):
    hull = ConvexHull(rnas[["x", "y"]])
    print(hull.area)
    return hull.area

area = data.groupby(level=0).aggregate(get_area).iloc[:, 0]
area.name = "area"

area.to_csv(snakemake.output[0], sep="\t", header=True)
