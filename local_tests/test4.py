from pertubations import UnitSphereDist
import numpy as np

samples = 1000

XYZ = np.zeros((samples, 3))


for jj in range(0, samples):
    XYZ[jj, :] = UnitSphereDist()


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(XYZ[:, 0], XYZ[:, 1], XYZ[:, 2])

plt.show()