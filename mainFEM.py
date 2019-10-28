#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 19:55:04 2019

@author: bjafib
"""

import numpy as np
from Q4Mesh import Q4Mesh
from Q4MeshPlot import PlotElemNumb, PlotNodeNumb, PlotDOFNumb
from Q4StiffnessMatrix import Q4GlobalStiffnessMatrix
from Displacement import Displacement
import matplotlib.pyplot as plt

# Order for Gauss quadrature
nGP = 2;

# Dimensions of beam
Length = 10
Height = 4
h = 1

# Number of elements in x and y direction
NElemX = 10
NElemY = 4

# Applied force
P = 0.001

# Constitutive matrix (for plane stress)
E = 1
nu = 0.3
Em = np.matrix([[1, nu, 0],
                 [nu, 1, 0],
                 [0, 0, (1-nu)/2]]) * E/(1-nu**2)

# Plot mesh
NElem, NNodes, NDOF, CoorX, CoorY, Nodes, DOFs = Q4Mesh(Length,Height,NElemX,NElemY)
#PlotElemNumb(CoorX,CoorY,NElem)
#PlotNodeNumb(CoorX,CoorY,Nodes,NNodes,NElem)
#PlotDOFNumb(CoorX,CoorY,DOFs,NDOF,NElem)

# Global stiffness matrix
K = Q4GlobalStiffnessMatrix(NDOF,NElem,CoorX,CoorY,DOFs,Em,h,nGP)

# Load vector
F = np.zeros((NDOF,1))
F[1,0] = -P # Node the applied load acting on

# Displacement
fixedDOFs = np.array([1,3,5,7,9,101,102])
freeDOFs = np.setdiff1d(DOFs, fixedDOFs).astype(int)
u, CoorXdef, CoorYdef = Displacement(F, K, fixedDOFs, freeDOFs, NDOF, CoorX, CoorY, DOFs)


#fig, ax = plt.subplots(1, 1)
#x, y = np.meshgrid(CoorXdef, CoorYdef)
#c = np.ones_like(x)
#ax.pcolormesh(x, y, c, facecolor='none', edgecolor='k')

x, y = np.meshgrid(CoorX, CoorY)
xdef, ydef = np.meshgrid(CoorXdef, CoorYdef)

#y, x = np.mgrid[:10, :10]
z = np.ones_like(x)

xdef, ydef = x**2, y**2 + x

fig, axes = plt.subplots(ncols=2)
axes[0].pcolormesh(x, y, z, facecolor='none', edgecolor='k')
axes[1].pcolormesh(xdef, ydef, z, facecolor='none', edgecolor='k')

axes[0].set(title='Original', xticks=[], yticks=[])
axes[1].set(title='Deformed', xticks=[], yticks=[])

plt.show()
