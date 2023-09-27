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
# print(bounding_distances)
# exit()


colors = cm.rainbow(np.linspace(0, 1, num_particles))
camera = Camera(plt.figure())
marker_sizes = 10*np.pi*radii**2

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
for frame in tqdm(range(num_frames)):
    plt.scatter(*trajectories[frame].T, c=colors, s=marker_sizes)
    for coords, bb in zip(trajectories[frame], bounding_distances):
        ax.add_patch(Rectangle(coords-bb*np.array([1.0, 1.0]), 2*bb, 2*bb, linewidth=1,
                               edgecolor="black", facecolor="none"))
    pairs = create_neighbor_pairs(neighbors_matrix[frame])
    for (i, j) in pairs:
        xi = trajectories[frame, i, 0]
        yi = trajectories[frame, i, 1]
        xj = trajectories[frame, j, 0]
        yj = trajectories[frame, j, 1]
        plt.plot([xi, xj], [yi, yj], "--", c="black", linewidth=1)
    camera.snap()

anim = camera.animate(blit=True)
FFwriter = animation.FFMpegWriter(
    fps=int(argv[4]),
    codec=None,
)
anim.save(f"{argv[2]}.mp4", writer=FFwriter)
