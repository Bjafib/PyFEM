#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 20:02:44 2019

@author: bjafib
"""

import numpy as np

def Q4Mesh(Length,Height,NElemX,NElemY):
    
    # Number of elements
    NElem = NElemX*NElemY
    # Total number of nodes
    NNodes=(NElemX+1)*(NElemY+1)
    # Total defrees of freedom
    NDOF = NNodes*2
    
    # Initiation of matrices
    CoorX = np.zeros((4,NElem))
    CoorY = np.zeros((4,NElem))
    Nodes = np.zeros((4,NElem))
    DOFs  = np.zeros((8,NElem))
    
    # Counting elements, nodes and dofs
    elem = 0 # Element counter
    for elemY in range(NElemY):
        # y-coordinates of element
        y1 = ((elemY)/NElemY)*Height
        y2 = y1
        y3 = ((elemY+1)/NElemY)*Height
        y4 = y3
        for elemX in range(NElemX):
            # x-coordinates of element
            x1 = ((elemX)/NElemX)*Length
            x2 = ((elemX+1)/NElemX)*Length
            x3 = x2
            x4 = x1
            # Element number
            elem = elem+1
            # Coordinates
            coorX = np.matrix([x1,x2,x3,x4])
            coorY = np.matrix([y1,y2,y3,y4])
            CoorX[0:4,elem-1:elem] = coorX.transpose()
            CoorY[0:4,elem-1:elem] = coorY.transpose()
            # Node numbers
            i1 = (NElemY+1)*(elemX)+elemY+1
            i2 = (NElemY+1)*(elemX+1)+elemY+1
            i3 = i2+1
            i4 = i1+1
            nodes = np.matrix([i1,i2,i3,i4])
            Nodes[0:4,elem-1:elem] = nodes.transpose()
            # DOF numbers
            dofs = np.matrix([(i1*2)-1, i1*2, (i2*2)-1, i2*2, (i3*2)-1, i3*2, (i4*2)-1, i4*2])
            DOFs[0:8,elem-1:elem] = dofs.transpose()
            
    return NElem, NNodes, NDOF, CoorX, CoorY, Nodes, DOFs