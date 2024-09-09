height = 20
width = 30
radius = 2

pds = PoissonDiskSampling(radius, height, width, seed=0)

import numpy as np
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

points = np.array(pds.points)
ax.scatter(points[:, 1], points[:, 2])
for p in points:
    ax.annotate(int(p[0]), (p[1], p[2]))
ax.set_aspect('equal')
plt.show()