import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
from matplotlib import animation
import numpy as np
from sys import argv
from celluloid import Camera
from tqdm import tqdm


data = np.load(argv[1])

# Process data
skip = int(argv[3])
trajectories = data["trajectories"][::skip, :, :]
num_particles = trajectories.shape[1]
num_frames = trajectories.shape[0]
width, height = data["box_size"]
neighbors_matrix = data["neighbors_matrix"]

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
check_neighbor_ids = [0, 10, 13, 27, 59]
focused_color = np.array([.0, .0, .0, 1.])
link_color = np.array([0., .0, 1., 1.])
for i in check_neighbor_ids:
    colors[i] = focused_color

for frame in tqdm(range(num_frames)):
    plt.xlim((0.0, width))
    plt.ylim((0.0, height))
    for idx in check_neighbor_ids:
        neighbor_ids = np.where(neighbors_matrix[frame, idx] > 0)[0]
        for n_id in neighbor_ids:
            p0x = trajectories[frame, idx, 0]
            p0y = trajectories[frame, idx, 1]
            p1x = trajectories[frame, n_id, 0]
            p1y = trajectories[frame, n_id, 1]
            plt.plot([p0x, p1x], [p0y, p1y],
                     c=link_color, linewidth=0.5)
    plt.scatter(*trajectories[frame].T, c=colors, s=25)
    camera.snap()
anim = camera.animate(blit=True)
FFwriter = animation.FFMpegWriter(
    fps=int(argv[4]),
    codec=None,
)
anim.save(f"{argv[2]}.mp4", writer=FFwriter)
