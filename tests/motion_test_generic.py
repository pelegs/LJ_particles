import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
from matplotlib import animation
import numpy as np
from sys import argv
from celluloid import Camera
from tqdm import tqdm


def create_neighbor_pairs(NM):
    Is, Js = np.where(NM == 1)
    pairs = set()
    for i, j in zip(Is, Js):
        pairs.add(frozenset({i, j}))
    return pairs


data = np.load(argv[1])

# Process data
skip = int(argv[3])
trajectories = data["trajectories"][::skip]
masses = data["masses"]
radii = data["radii"]
num_particles = trajectories.shape[1]
num_frames = trajectories.shape[0]
width, height = data["space_dimensions"]
neighbors_matrix = data["neighbors_matrix"][::skip]
bounding_distances = data["bounding_distances"]

if "AABB_min" in data:
    AABB_min = data["AABB_min"][::skip]
    AABB_max = data["AABB_max"][::skip]
    AABB_draw = True
else:
    AABB_draw = False

if "sort_by_x" in data:
    positions_sorted_x = data["sort_by_x"][::skip]
    positions_sorted_y = data["sort_by_y"][::skip]
    sort_particles = True
else:
    sort_particles = False
# print(positions_sorted_x.shape)
# print(positions_sorted_y.shape)
# exit()

colors = cm.rainbow(np.linspace(0, 1, num_particles))
camera = Camera(plt.figure())
marker_sizes = 10*np.pi*radii**2

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
for frame in tqdm(range(num_frames)):
    plt.scatter(*trajectories[frame].T, c=colors, s=marker_sizes)
    if AABB_draw:
        for pt_min, pt_max in zip(AABB_min[frame], AABB_max[frame]):
            bb = pt_max - pt_min
            ax.add_patch(Rectangle(pt_min, bb[0], bb[1], linewidth=1, linestyle="--",
                                   edgecolor="black", facecolor="none"))
    if sort_particles:
        id_min_x = positions_sorted_x[frame][0]
        id_min_y = positions_sorted_y[frame][0]
        pos_min_x = trajectories[frame, id_min_x, 0]
        pos_min_y = trajectories[frame, id_min_y, 1]
        plt.plot([pos_min_x, pos_min_x], [0, height],
                 color=colors[id_min_x])
        plt.plot([0, width], [pos_min_y, pos_min_y],
                 color=colors[id_min_y])
    camera.snap()

anim = camera.animate(blit=True)
FFwriter = animation.FFMpegWriter(
    fps=int(argv[4]),
    codec=None,
)
anim.save(f"{argv[2]}.mp4", writer=FFwriter)
