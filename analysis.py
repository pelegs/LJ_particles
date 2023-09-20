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
acum_disp = np.cumsum(displacments, axis=0)
sqr_acum_disp = acum_disp**2
mean_sqr_acum_disp = np.mean(sqr_acum_disp, axis=1)

# velocities = np.diff(trajectories, axis=0) / dt
# vels_norm = np.linalg.norm(velocities, axis=2)
# kinetic_energies = 0.5 * masses * vels_norm**2
# average_ke = np.mean(kinetic_energies, axis=1)

skip = 1000
plt.plot(times[::skip], displacments[::skip, :5])
plt.show()
