import matplotlib.pyplot as plt
# from matplotlib import cm
# from matplotlib.patches import Rectangle
# from matplotlib import animation
import numpy as np
from sys import argv


# Load data
data = np.load(argv[1])
dt = float(argv[2])

# Process data
trajectories = data["trajectories"][::, :, :]
masses = data["masses"]
radii = data["radii"]
num_particles = trajectories.shape[1]
num_steps = trajectories.shape[0]
width, height = data["box_size"]
times = np.arange(0.0, dt*num_steps, dt)

# Analysis
displacments = np.linalg.norm(trajectories, axis=2)
dis_diff = np.diff(displacments, axis=0)
diff_sqr = dis_diff**2
MSD = np.mean(diff_sqr, axis=1)

velocities = np.diff(trajectories, axis=0) / dt
vels_norm = np.linalg.norm(velocities, axis=2)
kinetic_energies = 0.5 * masses * vels_norm**2
total_ke = np.sum(kinetic_energies, axis=1)
average_ke = np.mean(kinetic_energies, axis=1)

skip = 1
# plt.plot(times[:-1:skip], total_ke[::skip])
plt.plot(times[:-1:skip], MSD[::skip])
plt.show()
