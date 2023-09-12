import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("LJ_test1.data")

r = data[300:600, 0]
U_LJ = data[300:600, 1]
F_LJ = data[300:600, 2]

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("Plot title")
ax1.set_xlabel("r")
ax1.set_ylabel("U(r)")

ax1.plot(r, U_LJ, c="r", label="Potential")
ax1.plot(r, F_LJ, c="b", label="Force")

leg = ax1.legend()

plt.show()
