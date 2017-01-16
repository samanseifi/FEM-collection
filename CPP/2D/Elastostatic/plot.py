#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

u = np.loadtxt('u.dat')
info = np.loadtxt('info.dat')

lx = info[0]
ly = info[1]
nx = info[2]
ny = info[3]

xx = np.linspace(0, lx, nx)
yy = np.linspace(0, ly, ny)
ux = np.reshape(u[0::2], (ny, nx))
uy = np.reshape(u[1::2], (ny, nx))

plt.contourf(xx, yy, ux, 100)
plt.show()
