#!/usr/bin/env python
from numpy import *
from matplotlib import pyplot as plt
from numpy.linalg import inv,det

import distmesh as dm
"""
 Uniform mesh class:
   Creates an ex by ey element mesh spanning a rectangular domain
   where 0 < x < lx and 0 < y < ly.
"""
class Mesh:
    # Initializes the mesh structure.
    def __init__(self, ex, ey, lx, ly):
        self.ex,self.ey = ex,ey
        self.lx,self.ly = lx,ly
        self.nx,self.ny = ex+1,ey+1
        self.hx,self.hy = lx/ex,ly/ey
        # Initialize nodal coordinates.
        self.nodes = []
        for y in linspace(0.0, ly, self.ny):
            for x in linspace(0.0, lx, self.nx):
                self.nodes.append([x,y])
        self.nodes = array(self.nodes)
        self.conn = []
        for j in range(self.ey):
            for i in range(self.ex):
                n0 = i + j*self.nx
                self.conn.append([n0, n0+1, n0+1+self.nx, n0+self.nx])
    # Returns the number of nodes.
    def num_nodes(self):    return self.nx*self.ny
    # Returns the number of elements.
    def num_elements(self): return self.ex*self.ey
    
class Hole:
	def __init__(self, lx, ly, cx, cy, r):
		fd = lambda p: dm.ddiff(dm.drectangle(p,-lx,lx,-ly,ly), dm.dcircle(p,cx,cx,r))
		fh = lambda p: 0.05+0.03*dm.dcircle(p,cx,cy,r)	
		p, t = dm.distmeshnd(fd, fh, 0.5, (-lx,-ly,lx,ly),[(-lx,-ly),(-lx,ly),(lx,-ly),(lx,ly)])
		self.nodes = p
		self.conn = t
	def num_nodes(self):
		size = self.nodes.shape
		return size[0]
	def num_elements(self):
		nel = self.conn.shape
		return nel[0]

def main(args):
    print 'constructing mesh...'
    mesh = Mesh(60,30,2.0,2.0)
    
    #mesh = Hole(20, 20, 0, 0, 5.0)
	

    # Plane-strain material tangent.
    E,v = 100.0,0.3
    C = E/(1.0+v)/(1.0-2.0*v) * array([[1.0-v, v,     0.0],
                                      [v,     1.0-v, 0.0],
                                      [0.0,   0.0,   0.5-v]])

    # Make stiffness matrix.
    K = zeros((2*mesh.num_nodes(), 2*mesh.num_nodes()))
    q4 = [[x/sqrt(3.0),y/sqrt(3.0)] for y in [-1.0,1.0] for x in [-1.0,1.0]]

    print 'assembling stiffness matrix...'
    B = zeros((3,8))
    for c in mesh.conn:
        xIe = mesh.nodes[c,:]
        Ke = zeros((8,8))
        for q in q4:
            dN = gradshape(q)
            J  = dot(dN, xIe).T
            dN = dot(inv(J), dN)
            # Assemble B matrix.
            B[0,0::2] = dN[0,:]
            B[1,1::2] = dN[1,:]
            B[2,0::2] = dN[1,:]
            B[2,1::2] = dN[0,:]

            Ke += dot(dot(B.T,C),B) * det(J)
        # Scatter operation.
        for i,I in enumerate(c):
            for j,J in enumerate(c):
                K[2*I,2*J]     += Ke[2*i,2*j]
                K[2*I+1,2*J]   += Ke[2*i+1,2*j]
                K[2*I+1,2*J+1] += Ke[2*i+1,2*j+1]
                K[2*I,2*J+1]   += Ke[2*i,2*j+1]
    # Assign nodal forces and boundary conditions.
    f = zeros((2*mesh.num_nodes()))
    for i in range(mesh.num_nodes()):
        if mesh.nodes[i,1] == 0.0:
            K[2*i,:]       = 0.0
            K[2*i+1,:]     = 0.0
            K[2*i,2*i]     = 1.0
            K[2*i+1,2*i+1] = 1.0
        if mesh.nodes[i,1] == mesh.ly:
            x = mesh.nodes[i,0]
            fbar = mesh.hx*(cos(8.0*pi*x/mesh.lx))
            fbar *= (x*(mesh.lx - x)) / mesh.lx**2
            f[2*i+1] = fbar
            if x == 0.0 or x == mesh.lx:
                f[2*i+1]  *= 0.5
    print 'solving linear system...'
    u = linalg.solve(K, f)
    print 'plotting displacement...'
    xx = linspace(0,mesh.lx,mesh.nx)
    yy = linspace(0,mesh.ly,mesh.ny)
    ux = reshape(u[0::2], (mesh.ny,mesh.nx))
    uy = reshape(u[1::2], (mesh.ny,mesh.nx))
    plt.contourf(xx, yy, uy**2+ux**2, 100)
    plt.show()

# Shape functions for a 4-node, isoparametric element.
def shape(xi):
    x,y = tuple(xi)
    N = [(1.0-x)*(1.0-y), (1.0+x)*(1.0-y), (1.0+x)*(1.0+y), (1.0-x)*(1.0+y)]
    return 0.25*array(N)

# Gradient of the shape functions for a 4-node, isoparametric element.
def gradshape(xi):
    x,y = tuple(xi)
    dN = [[-(1.0-y),  (1.0-y), (1.0+y), -(1.0+y)],
          [-(1.0-x), -(1.0+x), (1.0+x),  (1.0-x)]]
    return 0.25*array(dN)

if __name__ == '__main__':
    import sys
    main(sys.argv)

