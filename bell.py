
from matplotlib import pyplot as p
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from random import randint, uniform
from matplotlib import cm
import math
import random
import mpl_toolkits.mplot3d.art3d as art3d


def eucdist(a, b):
    assert len(a)==len(b), 'different dimensions'    
    return math.sqrt(sum((a[k] - b[k])**2 for k in range(len(a))))

#A REGION-BASED MULTI-ISSUE NEGOTIATION PROTOCOL FOR NONMONOTONIC UTILITY SPACES.pdf
def fbell(s, c, h, r):
    if eucdist(s, c) < r/2:     return h - 2 * h * (eucdist(s, c)**2) / (r**2)
    if eucdist(s, c) < r and eucdist(s, c) >=r/2:   return (2 * h / r ** 2) * (eucdist(s, c) - r )**2
    if eucdist(s, c) >= r:      return 0

# http://en.wikipedia.org/wiki/Gaussian_function
def f(x, y, h, o, c):
    a = 2
    return h * np.exp( - ((x- o[0])/(2*c**2)) ** 2 - ((y- o[1])/(2*c**2)) ** 2 )

fig = p.figure(figsize = (10, 8))
ax = fig.gca(projection='3d')
ax.set_aspect("auto")
ax.set_autoscale_on(True)
xlim, ylim, dx = -10., 10., .1

h = 2.
o = [0, 0, 0]
c = 1
folder = 'figures/'

ax.plot([o[0], o[0], ], [o[1], o[1]], [0, h], 'ro')

#width = 2 * c * np.sqrt(np.log(2))
#width = 2 * np.sqrt(2 * np.log(10)) * c
width = 2 * c * np.sqrt(np.log(2))

base = p.Circle((o[0], o[1]), width, facecolor='r', linewidth=0.6, alpha=.3)
ax.add_patch(base)
art3d.pathpatch_2d_to_3d(base, z=0, zdir='z')
    

number_of_points = 10
for point in range(number_of_points):
    pt = [uniform(-2, 2), uniform(-2, 2), uniform(0, h)]
    ax.plot([pt[0]], [pt[1]], [pt[2]], 'g^')
    print (' in cone     : ',  fbell(pt, o, h, width))
    print (' in function : ', f(pt[0], pt[1], h, o, c) >= pt[2])
    print ('\n')

###

x = y = np.arange(xlim, ylim, dx)
X, Y = np.meshgrid(x, y)
zs = np.array([ f(x, y, h, o, c) for x,y in zip(np.ravel(X), np.ravel(Y))])

Z = zs.reshape(X.shape)
surf = ax.plot_surface(X, Y, Z, alpha=0.02, linewidth=0.15, edgecolors='b') #, cmap=cm.jet)
#fig.colorbar(surf, shrink=0.5, aspect=5)


ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("f")
ax.set_xlim([xlim, ylim])
ax.set_ylim([xlim, ylim])
ax.set_zlim([0., 2.5])

p.savefig(folder + 'Bell.pdf', format='pdf', dpi=1000)

