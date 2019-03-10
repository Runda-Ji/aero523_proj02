# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 16:45:25 2018

@author: rundaji
"""
import numpy as np
import matplotlib.pyplot as plt
from flux import HLLE_flux, flux

def vary_norm_vec(M_L,M_R,title):
    #set two states as constants, check the impact of the norm vector
    global gamma, alpha_L, alpha_R;
    U_L = np.array([1,
                    M_L*np.cos(alpha_L),
                    M_L*np.sin(alpha_L),
                    1/((gamma-1)*gamma)+0.5*M_L**2]);
    U_R = np.array([1,
                    M_R*np.cos(alpha_R),
                    M_R*np.sin(alpha_R),
                    1/((gamma-1)*gamma)+0.5+M_R**2]);
    #--------------------------------------------------------------------------
    N = 1000;
    theta = [None]*N;
    F = [None]*N;
    for i in range(0, N):
        theta[i] =0.5* i/N*np.pi; #theta in [0,pi/2]
        n = np.array([np.cos(theta[i]),np.sin(theta[i])]);
        theta[i] = theta[i]/np.pi*180;
        F[i],s = HLLE_flux(U_L,U_R,n);
    #--------------------------------------------------------------------------
    f = plt.figure(figsize=(8,6));
    l1,l2,l3,l4 = plt.plot(theta,F,'-');    
    plt.xlabel(r'Normal vector direction, ${\theta}$ (degree)',fontsize = 16);
    plt.ylabel('Flux',fontsize = 16);
    plt.grid();
    plt.legend([l1,l2,l3,l4], ["Mass flux","Momentum flux, x-component","Momentum flux, y-component","Energy flux"], loc=1);
    plt.savefig('..\\figure_prob_1\\vary_norm_vec_%s.pdf' %title, dpi=150);    
    plt.show();
    plt.close(f);
    return 0;

def vary_right_Mach_number(M_L,title):
    #set the left Mach # as constant, check the impact of changing right Mach #
    global gamma, alpha_L, alpha_R;
    n = np.array([1,0]);
    U_L = np.array([1,
                    M_L*np.cos(alpha_L),
                    M_L*np.sin(alpha_L),
                    1/((gamma-1)*gamma)+0.5*M_L**2]);
    #--------------------------------------------------------------------------
    N = 30;
    M_R = [None]*N
    F = [None]*N;
    for i in range(0,N):
        M_R[i] = (i+1)/10; # M_R in [0.1,3.1]
        U_R = np.array([1,
                        M_R[i]*np.cos(alpha_R),
                        M_R[i]*np.sin(alpha_R),
                        1/((gamma-1)*gamma)+0.5*M_R[i]**2]);
        F[i],s = HLLE_flux(U_L,U_R,n);    
    #--------------------------------------------------------------------------
    f = plt.figure(figsize=(8,6));
    l1,l2,l3,l4 = plt.plot(M_R,F,'-');
    plt.xlabel(r'$M_R$',fontsize = 16);
    plt.ylabel('Flux',fontsize = 16);
    plt.grid();
    plt.legend([l1,l2,l3,l4], ["Mass flux","Momentum flux, x-component","Momentum flux, y-component","Energy flux"], loc=0);
    plt.savefig('..\\figure_prob_1\\vary_right_Mach_number_%s.pdf' %title, dpi=150);
    plt.show();
    plt.close(f);
    return 0;

def compare_with_Euler():
    global gamma;
    alpha = 0;
    n = np.array([1,0]);
    #--------------------------------------------------------------------------
    N = 30;
    M = [None]*N;
    F_Euler = [None]*N;
    F_HLLE = [None]*N;
    #--------------------------------------------------------------------------
    for i in range(0,N):
        M[i] = (i+1)/10; # M_R in [0.1,3.1]
        U = np.array([1,
                      M[i]*np.cos(alpha),
                      M[i]*np.sin(alpha),
                      1/((gamma-1)*gamma)+0.5*M[i]**2]);
        prop = flux(U);
        F_Euler[i] = prop['F'][0]*n[0] + prop['F'][1]*n[1];
        F_HLLE[i],s = HLLE_flux(U,U,n);
    #--------------------------------------------------------------------------
    f = plt.figure(figsize=(8,6));
    l1,l2,l3,l4 = plt.plot(M,F_Euler,'x');
    l5,l6,l7,l8 = plt.plot(M,F_HLLE,'-');
    plt.xlabel(r'$M_L=M_R=M$',fontsize = 16);
    plt.ylabel('Flux',fontsize = 16);
    plt.grid();
    plt.legend([l1,l2,l3,l4,l5,l6,l7,l8], ["Mass flux, Euler equation", \
                                           "Momentum flux, x-component, Euler equation", \
                                           "Momentum flux, y-component, Euler equation", \
                                           "Energy flux, Euler equation", \
                                           "Mass flux, HLLE ", \
                                           "Momentum flux, x-component, HLLE", \
                                           "Momentum flux, y-component, HLLE", \
                                           "Energy flux, HLLE"], loc=0);
    plt.savefig('..\\figure_prob_1\\compare_with_Euler.pdf',dpi=150);
    plt.show();
    plt.close(f);
    return 0;

def main():
    global gamma, alpha_L, alpha_R;
    gamma = 1.3; alpha_L = 0; alpha_R = 0;
    ###########################################################################
    #subsonic case
    M_L = 0.06; M_R = 0.05;
    vary_norm_vec(M_L,M_R,'Subsonic');
    #supersonic case
    M_L = 2.06; M_R = 2.05
    vary_norm_vec(M_L,M_R,'Supersonic');
    ###########################################################################
    #subsonic case
    M_L = 0.5;
    vary_right_Mach_number(M_L,'Subsonic');
    #supersonic case
    M_L = 1.3;
    vary_right_Mach_number(M_L,'Supersonic');
    ###########################################################################
    compare_with_Euler();
    return 0;
    


if __name__=="__main__":
    main()