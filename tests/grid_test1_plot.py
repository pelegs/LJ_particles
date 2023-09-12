import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


data = np.loadtxt("tests/grid_test1.data", comments="#")
with open("tests/grid_test1.data", "r") as f:
    params = [next(f)
              for _ in range(1)][0].strip("# ").split(" ")
w, h = [float(x) for x in params[:2]]
nx, ny = [int(x) for x in params[2:]]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title("Motion test")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xticks(np.linspace(0, w, nx+1))
ax.set_yticks(np.linspace(0, h, ny+1))
ax.scatter(data[:, 0], data[:, 1])

txt_offset = np.array([5.0, 0.0])
for i, txt in enumerate(data[:, 2:].astype(int)):
    ax.annotate(txt, data[i, :2]+txt_offset)
ax.add_patch(Rectangle((0.0, 0.0), w, h, linewidth=4,
                       edgecolor='b', facecolor='none'))

plt.grid()
plt.show()
