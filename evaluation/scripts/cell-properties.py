import pandas as pd
from scipy.spatial import ConvexHull

data = pd.read_table(snakemake.input[0], index_col=0)

def get_area(rnas):
    hull = ConvexHull(rnas[["x", "y"]])
    return hull.area

def get_first(rnas):
    return rnas.iloc[0]

cells = data.groupby(level=0)

area = cells[["x", "y"]].aggregate(get_area).iloc[:, 0]
x = cells["cell_x"].aggregate(get_first)
y = cells["cell_y"].aggregate(get_first)

props = pd.DataFrame({"area": area, "x": x, "y": y})

props.to_csv(snakemake.output[0], sep="\t")
