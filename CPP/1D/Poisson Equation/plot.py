#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np


# FEM Solution
f_fem = np.loadtxt('f.txt')
m = np.size(f_fem)
x1 = np.linspace(0, 1, m)

# Exact Solution
x2 = np.linspace(0, 1, 1000)
f_exact = - (x2**4)/12.0 + (x2**3)/3.0 - x2/4.0

# Plotting
plt.plot(x1, f_fem, '-o', label='Exact')
plt.plot(x2, f_exact, label='FEM')

plt.legend(loc='upper left', fontsize = 12)
plt.xlabel('$x$', fontsize = 14)
plt.ylabel('$u(x)$', fontsize = 14)
plt.savefig('f.pdf')

plt.show()
