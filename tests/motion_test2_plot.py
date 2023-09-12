import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from celluloid import Camera
from tqdm import tqdm


data = np.loadtxt("tests/motion_test2.data", comments="#")
num_particles = data[0].shape[0]
num_steps = data.shape[0]//2
data = data.reshape((num_steps, 2, num_particles))

# numpoints = 10
# points = np.random.random((2, numpoints))
colors = cm.rainbow(np.linspace(0, 1, num_particles))
camera = Camera(plt.figure())
for step in tqdm(range(num_steps)):
    plt.scatter(*data[step], c=colors, s=400)
    camera.snap()
anim = camera.animate(blit=True)
anim.save("tests/motion_test2.mp4", fps=60)
