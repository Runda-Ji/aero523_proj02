# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 11:35:28 2018

@author: rundaji
"""

import numpy as np

def initialization(mesh,alpha):
    gamma = 1.3;
    M_0 = 1.3;
    U_0 = np.array([1,
                    M_0*np.cos(alpha),
                    M_0*np.sin(alpha),
                    1/((gamma-1)*gamma)+0.5*M_0**2]);
    
    for i in range(0,mesh['nElem']):
        mesh['Elems'][i].state = U_0;
    return 0;

def read_state(mesh,fname):
    f = open(fname, "r");
    for i in range(0,mesh['nElem']):
        mesh['Elems'][i].state = [float(string) for string in f.readline().split()];
        mesh['Elems'][i].state = np.asarray(mesh['Elems'][i].state);
    f.close();
    return 0;

def find_U_inf(alpha):
    M_inf = 8;
    alpha = alpha/180*np.pi;
    gamma = 1.3;
    U_inf = np.array([1,
                      M_inf*np.cos(alpha),
                      M_inf*np.sin(alpha),
                      1/((gamma-1)*gamma)+0.5*M_inf**2]);
    return U_inf;
