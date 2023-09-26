import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from sys import argv


data = np.load(argv[1])

# Process data
trajectories = data["trajectories"][::, :, :]
num_particles = trajectories.shape[1]
num_frames = trajectories.shape[0]
width, height = data["space_dimensions"]
sorted_by_x = data["sort_by_x"]
sorted_by_y = data["sort_by_y"]
print(sorted_by_x)
print(sorted_by_y)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title("Motion test")
ax.set_xlabel("x")
ax.set_ylabel("y")
# ax.set_xticks(np.linspace(0, width, 10))
# ax.set_yticks(np.linspace(0, height, 10))
ax.add_patch(Rectangle((0.0, 0.0), width, height, linewidth=4,
                       edgecolor='b', facecolor='none'))
ax.scatter(*trajectories[0].T)
ds = 10
sq = np.array([ds, ds])
for i, pos in enumerate(trajectories[0]):
    ax.add_patch(Rectangle(pos-sq, 2*ds, 2*ds, linewidth=2,
                 edgecolor="r", facecolor="none"))
    ax.annotate(f"{i}", pos, fontsize=25)

plt.show()
