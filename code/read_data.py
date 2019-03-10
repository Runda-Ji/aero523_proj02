# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 21:31:12 2018

@author: Runda Ji
"""

from my_class import CELL, BOUNDARY
import numpy as np

def read_data(fname):
    f = open(fname, "r");
    #read some general info
    nNode,nElem,Dim = [int(string) for string in f.readline().split()];
    node_pos = [None]*nNode;
    #read the position of nodes
    for i in range(0,nNode):
        x,y = [float(string) for string in f.readline().split()];
        node_pos[i] = [x,y];
    node_pos =  np.asarray(node_pos);
    #read the number of boundary groups
    nBGroup = int(f.readline());
    #read the boundaries
    boundary = [None]*nBGroup;
    for i in range(0,nBGroup):
        nBFace,nf,Title = f.readline().split();
        nBFace = int(nBFace);
        nf = int(nf);
        B = [None]*nBFace;
        for j in range(0,nBFace):
            n_0, n_1 = f.readline().split();
            # the given index start from 1, we want the index start from 0
            n_0 = int(n_0)-1;
            n_1 = int(n_1)-1;
            if n_0 > n_1:
                n_0, n_1 = n_1, n_0;
            B[j] = {'n0': n_0, 'n1': n_1};
        boundary[i] = BOUNDARY(nBFace, nf, Title, B);
    #read cell info
    nElem,Order,Basis = f.readline().split();
    nElem = int(nElem);
    Order = int(Order);
    E = [None]*nElem;
    for i in range(0,nElem):
        v_0,v_1,v_2 = [int(string) for string in f.readline().split()];
        # the given index start from 1, we want the index start from 0
        v_0 = v_0-1;
        v_1 = v_1-1;
        v_2 = v_2-1;
        e_2 = {'n0': v_0, 'n1': v_1, 't': i, 'e':2};
        e_0 = {'n0': v_1, 'n1': v_2, 't': i, 'e':0};
        e_1 = {'n0': v_2, 'n1': v_0, 't': i, 'e':1};
        vertex = [v_0, v_1, v_2];
        edges = [e_0, e_1, e_2];
        # vertex, edge, tri, A, adj_cell, state, R, dt
        E[i] = CELL(vertex, edges, i, [], None, None, None);
    mesh = {'nNode':nNode, 'nElem':nElem, 'node_pos':node_pos, 'nBGroup':nBGroup, 'boundary':boundary, 'Elems':E, 'nEdge': 0, 'Edges':[]};
    f.close();
    return mesh;