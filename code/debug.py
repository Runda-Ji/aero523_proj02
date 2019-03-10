# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt


x = mesh['node_pos'][:,0];
y = mesh['node_pos'][:,1];
tri = [None]*mesh['nElem'];
for i in range(0,mesh['nElem']):
    tri[i] = mesh['Elems'][i].vertex;
f2 = plt.figure(figsize=(8,6));
plt.triplot(x, y, tri, 'k-', lw=0.5);
pos = [];
pos.append([0.833058,0.119717]);
pos.append([0.800000,0.100000]);
pos.append([0.800000,0.150000]);
pos = np.asarray(pos);
plt.scatter(pos[:,0],pos[:,1],s=5,c='r');
plt.axis('equal');
plt.xlim((-1, 1));
plt.ylim((-1, 1));