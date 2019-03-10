# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 21:13:12 2018

@author: Runda Ji
"""

class EDGE():
    def __init__(self, node_0, node_1, t1, e1, t2, e2, norm_vec, length, s):
        self.node_0 = node_0;
        self.node_1 = node_1;
        self.t1 = t1;
        self.e1 = e1;
        self.t2 = t2;
        self.e2 = e2;
        self.norm_vec = norm_vec;
        self.length = length;
        self.s = s;

class CELL():
    def __init__(self, vertex, edge, tri, adj_cell, state, R, dt):
        self.vertex = vertex;
        self.edge = edge;
        self.tri = tri; #the No. of triangle (cell)
        self.adj_cell = adj_cell;
        self.state = state;
        self.R = R;
        self.dt = dt;
        
        
class BOUNDARY():
    def __init__(self, nBFace, nf, Title, B):
            self.nBFace = nBFace;
            self.nf = nf;
            self.Title = Title;
            self.B = B;

