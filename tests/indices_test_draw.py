import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np

with open("tests/indices_test1.txt", "r") as f:
    data = f.readlines()

nx, ny = [int(x) for x in data[0].split(" ")]
cell_i, cell_j = [int(x) for x in data[1].split(" ")]
data = np.loadtxt("tests/indices_test1.txt", skiprows=2).astype(int)
mat = np.zeros((nx, ny))
mat[cell_i, cell_j] = 2.0
for (i, j) in data:
    mat[i, j] = 1.0

# create discrete colormap
cmap = colors.ListedColormap(["red", "blue", "green"])
bounds = [0, 1, 2]
norm = colors.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots()
ax.imshow(mat, cmap=cmap)

# draw gridlines
ax.grid(which="major", axis="both", linestyle="-", color="k", linewidth=2)
ax.set_xticks(np.arange(-0.5, nx, 1))
ax.set_yticks(np.arange(-0.5, ny, 1))

plt.show()
