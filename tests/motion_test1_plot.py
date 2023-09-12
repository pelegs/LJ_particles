import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("tests/motion_test1.data")

r = data[:, 0]
p1 = data[:, 1]
p2 = data[:, 2]

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("Motion test")
ax1.set_xlabel("t")
ax1.set_ylabel("x")

ax1.plot(r, p1, c="r", label="p1 pos")
ax1.plot(r, p2, c="b", label="p2 pos")

leg = ax1.legend()

plt.show()
