import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
import colorsys


def rand_cmap(nlabels, type="bright", first_color_black=True, last_color_black=False, verbose=True):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    """
    if type not in ("bright", "soft"):
        print("Please choose 'bright' or 'soft' for type")
        return

    if verbose:
        print("Number of labels: " + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == "bright":
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                          np.random.uniform(low=0.2, high=1),
                          np.random.uniform(low=0.9, high=1)) for i in range(nlabels)]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(
                HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list(
            "new_map", randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == "soft":
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high)) for i in range(nlabels)]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list(
            "new_map", randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing="proportional", ticks=None,
                                   boundaries=bounds, format="%1i", orientation=u"horizontal")

    return random_colormap


data = np.loadtxt("tests/indices_test3.data", comments="#")
with open("tests/indices_test3.data", "r") as f:
    params = f.readline().strip("#  \n").split(" ")
    neighbors = f.readline().split(": ")
n_rows, n_cols = [int(x) for x in params[:2]]
width, height = [float(x) for x in params[2:4]]
neighbor_ids = [int(x) for x in neighbors[-1].strip().split(" ")]
num_particles = data.shape[1]

new_cmap = rand_cmap(n_rows*n_cols, type="bright",
                     first_color_black=False, last_color_black=False, verbose=False)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title("Indices test")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xticks(np.linspace(0, width, n_cols+1))
ax.set_yticks(np.linspace(0, height, n_rows+1))
ax.set_aspect('equal', adjustable='box')

for id in neighbor_ids:
    # print(data[:2, id])
    xs = np.array([data[0, id], data[0, 137]])
    ys = np.array([data[1, id], data[1, 137]])
    ax.plot(xs, ys, '-', c="grey")

ax.scatter(data[0], data[1], c=data[2], cmap=new_cmap)

# txt_offset = np.array([-1.0, 0.5])
# for i, txt in enumerate(data[2].astype(int)):
#     ax.annotate(txt, data[:2, i]+txt_offset)
# ax.add_patch(Rectangle((0.0, 0.0), width, height, linewidth=4,
#                        edgecolor="b", facecolor="none"))

plt.grid()
plt.show()
