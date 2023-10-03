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
# neighbors_matrix = data["neighbors_matrix"][::skip]
# bounding_distances = data["bounding_distances"]
# forces = data["forces"]
walls = data["walls_data"]

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

AABB_draw = sort_particles = annotate_sort = False
draw_extreme_lines = draw_neighbors = draw_force = True

colors = cm.rainbow(np.linspace(0, 1, num_particles))
camera = Camera(plt.figure())
marker_sizes = 0.25*np.pi*radii**2

# plt.figsize = (1000*px, 1000*px)
plt.title("Neighbors test using sort-n-sweep")
plt.xlabel = "x"
plt.ylabel = "y"
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
lbl_dx = np.array([3.0, 0.0])
lbl_dy = np.array([0.0, 5.0])
for frame in tqdm(range(num_frames)):
    plt.xlim(0.0, width)
    plt.ylim(0.0, height)
    plt.scatter(*trajectories[frame].T, c=colors, s=marker_sizes)
    for wall in walls:
        plt.plot(*wall.T, color="black", linewidth=3)
    # if draw_force:
    #     pass
    #
    # if AABB_draw:
    #     for pt_min, pt_max in zip(AABB_min[frame], AABB_max[frame]):
    #         bb = pt_max - pt_min
    #         ax.add_patch(Rectangle(pt_min, bb[0], bb[1], linewidth=0.5, linestyle="--",
    #                                edgecolor="black", facecolor="none"))
    # if sort_particles:
    #     if draw_extreme_lines:
    #         id_min_x = positions_sorted_x[frame][0]
    #         id_min_y = positions_sorted_y[frame][0]
    #         pos_min_x = trajectories[frame, id_min_x,
    #                                  0] - bounding_distances[id_min_x]
    #         pos_min_y = trajectories[frame, id_min_y,
    #                                  1] - bounding_distances[id_min_y]
    #         plt.plot([pos_min_x, pos_min_x], [0, height],
    #                  color=colors[id_min_x])
    #         plt.plot([0, width], [pos_min_y, pos_min_y],
    #                  color=colors[id_min_y])
    #
    #         id_max_x = positions_sorted_x[frame][-1]
    #         id_max_y = positions_sorted_y[frame][-1]
    #         pos_max_x = trajectories[frame, id_max_x,
    #                                  0] + bounding_distances[id_max_x]
    #         pos_max_y = trajectories[frame, id_max_y,
    #                                  1] + bounding_distances[id_max_y]
    #         plt.plot([pos_max_x, pos_max_x], [0, height],
    #                  color=colors[id_max_x])
    #         plt.plot([0, width], [pos_max_y, pos_max_y],
    #                  color=colors[id_max_y])
    #
    #     for i, id_x in enumerate(positions_sorted_x[frame]):
    #         pos = trajectories[frame, id_x]
    #         if annotate_sort:
    #             plt.annotate(
    #                 f"({i},", pos-lbl_dx+lbl_dy)
    #         if draw_neighbors:
    #             neighbor_ids = np.where(neighbors_matrix[frame][id_x] == 1)[0]
    #             for n_id in neighbor_ids:
    #                 pos_neighbor = trajectories[frame, n_id]
    #                 x_values = [pos[0], pos_neighbor[0]]
    #                 y_values = [pos[1], pos_neighbor[1]]
    #                 plt.plot(x_values, y_values,
    #                          linestyle="--", color=colors[id_x])
    #             # neighbors_list_txt = ",".join(
    #             #     map(str, neighbor_ids))
    #             # plt.annotate(f"{id_x}: {neighbors_list_txt}", pos)
    #     if annotate_sort:
    #         for j, id_y in enumerate(positions_sorted_y[frame]):
    #             pos = trajectories[frame, id_y]
    #             plt.annotate(f"{j})",
    #                          pos+lbl_dx+lbl_dy)
    camera.snap()

anim = camera.animate(blit=True)
FFwriter = animation.FFMpegWriter(
    fps=int(argv[4]),
    codec=None,
)
anim.save(f"{argv[2]}.mp4", writer=FFwriter)
