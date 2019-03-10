# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 16:43:53 2018

@author: rundaji
"""
import numpy as np

def flux(U):
    F_x = [None]*4;
    F_y = [None]*4;
    gamma = 1.3;
    rho = U[0];
    u = U[1]/U[0];
    v = U[2]/U[0];
    E = U[3]/U[0];
    p = (gamma-1)*(rho*E-0.5*rho*(u**2+v**2));
    #--------------------------------------------------------------------------
    if p < 0:
        p = np.abs(p);
        print('Error! Interior negative pressure!');
        print('Take the absolute value of pressure!');
        #return 'Negative Pressure';
    #--------------------------------------------------------------------------
    a = np.sqrt(gamma*p/rho);
    H = E + p/rho;
    
    F_x[0] = rho*u;
    F_x[1] = rho*u**2+p;
    F_x[2] = rho*u*v;
    F_x[3] = rho*u*H;
    F_y[0] = rho*v;
    F_y[1] = rho*v*u;
    F_y[2] = rho*v**2+p;
    F_y[3] = rho*v*H;
    
    V = np.array([u,v]);
    F = np.array([F_x,F_y]);
    #key properties of the state
    prop = {'V': V, 'F': F, 'a': a};
    return prop;
    
def HLLE_flux(U_L,U_R,n):
    #U is state
    #V = (u,v) is velocity
    prop_L = flux(U_L);
    prop_R = flux(U_R);
    
    #--------------------------------------------------------------------------
    if prop_L == 'Negative Pressure' or prop_R == 'Negative Pressure':
        print('Unable to compute HLLE flux');
        return -1,-1;
    #--------------------------------------------------------------------------
    
    Vn_L = np.dot(prop_L['V'],n);
    Vn_R = np.dot(prop_R['V'],n);
    F_L = prop_L['F'][0]*n[0] + prop_L['F'][1]*n[1];
    F_R = prop_R['F'][0]*n[0] + prop_R['F'][1]*n[1];
    c_L = prop_L['a'];
    c_R = prop_R['a'];
    
    s_L_min = np.minimum(0,Vn_L-c_L);
    s_L_max = np.maximum(0,Vn_L+c_L);
    s_R_min = np.minimum(0,Vn_R-c_R);
    s_R_max = np.maximum(0,Vn_R+c_R);
    s_min = np.minimum(s_L_min,s_R_min);
    s_max = np.maximum(s_L_max,s_R_max);
    s = np.maximum(np.abs(Vn_L)+c_L,np.abs(Vn_R)+c_R); # for computing time step
    F = 0.5*(F_L+F_R) \
       -0.5*(s_max+s_min)/(s_max-s_min)*(F_R-F_L) \
       +(s_max*s_min)/(s_max-s_min)*(U_R-U_L);
    return F,s;

def wall_flux(U,n):
    gamma = 1.3;
    rho = U[0];
    u = U[1]/U[0];
    v = U[2]/U[0];
    E = U[3]/U[0];
    V = np.array([u,v]);
    p = (gamma-1)*(rho*E-0.5*rho*(V[0]**2+V[1]**2));
    V_b = V - np.dot(V,n)*n;
    p_b = (gamma-1)*(rho*E-0.5*rho*(V_b[0]**2+V_b[1]**2));
    F_b = np.array([0,\
                    p_b*n[0],\
                    p_b*n[1],\
                    0]);
    
    #--------------------------------------------------------------------------
    if p < 0:
        p = np.abs(p);
        print('Error! Negative pressure on the wall!');
        print('Take the absolute value of pressure!');
        #return -1,-1;
    #--------------------------------------------------------------------------
    
    s = np.sqrt(gamma*p/rho);
    return F_b,s;