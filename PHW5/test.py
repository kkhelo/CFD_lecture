
import numpy as np 
import matplotlib.pyplot as plt
P = np.load('P.npy')

max = np.max(P)
min = np.min(P)
step = (max - min) / 50

m = 20
dx = dy = 1/m
x = np.arange(0 - 1*dx, 1 + 2*dx, dx)
y = np.arange(1 + 1*dy, 0 - 2*dy, -dy)
X, Y = np.meshgrid(x, y)

plt.figure()
contour = plt.contourf(X, Y, P, levels = np.arange(min, max + step, step))
plt.colorbar(contour)
plt.savefig('test.png')
plt.show()
