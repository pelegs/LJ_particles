import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
import numpy as np
from celluloid import Camera
from tqdm import tqdm


data = np.loadtxt("tests/motion_test3.data", comments="#")
num_particles = data[0].shape[0]
num_steps = data.shape[0]//2
data = data.reshape((num_steps, 2, num_particles))

width, height = 700.0, 700.0

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
for step in tqdm(range(num_steps)):
    plt.xlim((0.0, width))
    plt.ylim((0.0, height))
    plt.scatter(*data[step], c=colors, s=100)
    camera.snap()
anim = camera.animate(blit=True)
anim.save("tests/motion_test3.mp4", fps=60)
