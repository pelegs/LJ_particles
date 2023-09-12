import numpy as np
from scipy.stats import maxwell
import matplotlib.pyplot as plt


def maxw(size=None):
    """Generates size samples of maxwell"""
    vx = np.random.normal(size=size)
    vy = np.random.normal(size=size)
    vz = np.random.normal(size=size)
    return np.sqrt(vx*vx + vy*vy + vz*vz)


if __name__ == "__main__":
    # mdata = maxw(100000)
    # h, bins = np.histogram(mdata, bins=101, range=(0.0, 10.0))
    #
    # x = np.linspace(0.0, 10.0, 100)
    # rv = maxwell()
    #
    num_particles = 1000
    vels_q1 = np.random.uniform(size=(num_particles, 2))
    for i, v in enumerate(vels_q1):
        vels_q1[i] = v/np.linalg.norm(v)
    vels_q2 = vels_q1.copy()
    vels_q2[:, 0] *= -1.0
    vels_q3 = vels_q2.copy()
    vels_q3[:, 1] *= -1.0
    vels_q4 = vels_q3.copy()
    vels_q4[:, 0] *= -1.0
    vels = np.concatenate((vels_q1, vels_q2, vels_q3, vels_q4), axis=0)
    vel_cm = np.mean(vels, axis=0)
    print(vel_cm, np.linalg.norm(vel_cm))
    # origin = np.zeros(2)
    #
    # fig, ax = plt.subplots(1, 1)
    # # ax.hist(mdata, bins=bins, density=True)
    # ax.scatter(vels[:, 0], vels[:, 1], color="red", label="whatever")
    # # plt.quiver(*origin, vel_cm[0], vel_cm[1], color=[1, 0, 0])
    # plt.title("Something")
    # plt.show()
