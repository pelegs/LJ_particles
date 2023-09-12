import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("tests/indices_test3.data", comments="#")
with open("tests/indices_test3.data", "r") as f:
    params = f.readline().strip("#  \n").split(" ")
n_rows, n_cols = [int(x) for x in params[:2]]
width, height = [float(x) for x in params[2:]]
num_particles = data.shape[1]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title("Indices test")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xticks(np.linspace(0, width, n_cols+1))
ax.set_yticks(np.linspace(0, height, n_rows+1))
ax.scatter(data[0], data[1], c=data[2], cmap=new_cmap)

# txt_offset = np.array([-1.0, 0.5])
# for i, txt in enumerate(data[2].astype(int)):
#     ax.annotate(txt, data[:2, i]+txt_offset)
# ax.add_patch(Rectangle((0.0, 0.0), width, height, linewidth=4,
#                        edgecolor="b", facecolor="none"))

plt.grid()
plt.show()
