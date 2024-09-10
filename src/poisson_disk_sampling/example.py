import matplotlib.pyplot as plt
import numpy as np

import poisson_disk_sampling.PoissonDiskSampling as pds

height = 20
width = 30
radius = 2

samp = pds.sample_points(radius, width, height, seed=0)
points = np.array(samp)

fig, ax = plt.subplots()
ax.scatter(points[:, 1], points[:, 2])
for p in points:
    ax.annotate(int(p[0]), (p[1], p[2]))
ax.set_aspect('equal')
plt.show()
