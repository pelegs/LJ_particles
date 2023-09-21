import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
from matplotlib import animation
import numpy as np
from sys import argv
from celluloid import Camera
from tqdm import tqdm
from itertools import combinations


data = np.load(argv[1])

# Process data
skip = int(argv[3])
trajectories = data["trajectories"][::skip, :, :]
# masses = data["masses"]
# radii = data["radii"]
num_particles = trajectories.shape[1]
num_frames = trajectories.shape[0]
width, height = data["box_size"]
# neighbors_matrix = data["neighbors_matrix"][::skip, :, :]
# for i, frame in enumerate(neighbors_matrix):
#     if not np.array_equal(frame, frame.T):
#         print(i)
# exit()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title("Motion test")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xticks(np.linspace(0, width, 10))
ax.set_yticks(np.linspace(0, height, 10))
ax.add_patch(Rectangle((0.0, 0.0), width, height, linewidth=4,
                       edgecolor='b', facecolor='none'))

colors = cm.rainbow(np.linspace(0, 1, num_particles))
camera = Camera(plt.figure())
connected_ids = [20, 21, 28, 29]
connected_color = np.array([.0, .0, .0, 1.])
link_color = np.array([0., .0, 1., 1.])
for i in connected_ids:
    colors[i] = connected_color

# set marker sizes
# marker_sizes = radii*50
marker_sizes = 50

for frame in tqdm(range(num_frames)):
    plt.xlim((0.0, width))
    plt.ylim((0.0, height))
    for (id1, id2) in combinations(connected_ids, 2):
        xs = np.array([trajectories[frame, id1, 0],
                      trajectories[frame, id2, 0]])
        ys = np.array([trajectories[frame, id1, 1],
                      trajectories[frame, id2, 1]])
        plt.plot(xs, ys, "--", connected_color, linewidth=1)
    plt.scatter(*trajectories[frame].T, c=colors, s=marker_sizes)
    camera.snap()
anim = camera.animate(blit=True)
FFwriter = animation.FFMpegWriter(
    fps=int(argv[4]),
    codec=None,
)
anim.save(f"{argv[2]}.mp4", writer=FFwriter)
