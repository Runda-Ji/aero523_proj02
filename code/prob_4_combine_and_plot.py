# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 10:23:21 2018

@author: rundaji
"""

#from my_class import EDGE, CELL, BOUNDARY

from read_data import read_data
from pre_calculation import pre_calculation
from initialization import find_U_inf,read_state
from time_iteration import find_coefs
from post_processing import plot_mesh, plot_pressure_coef, plot_mach

import matplotlib.pyplot as plt

def main():
    global mesh,coefs,nElem,coef_p;
    alpha = 5;
    U_inf = find_U_inf(alpha);
    coefs = [];
    nElem = [];
    coef_p = [];
    folder = ('figure_prob_4\\alpha_%d' %alpha);
    
    mesh = read_data('..\\mesh\\capsule.gri');
    nElem.append(mesh['nElem']);
    pre_calculation(mesh);
    read_state(mesh,'..\\state_after_iteration\\Base_line.txt');
    cl = find_coefs(mesh,U_inf)['cl'];
    cd = find_coefs(mesh,U_inf)['cd'];
    cm = find_coefs(mesh,U_inf)['cm'];
    coefs.append([cl,cd,cm]);
    p_versus_theta = plot_pressure_coef(mesh,folder,'Base_line_cp_versus_theta',U_inf);
    coef_p.append(p_versus_theta);
    
    
    for i in range(1,6):
        mesh = read_data('..\\mesh\\alpha_%d_ada_iter_%d.gri' %(alpha,i));
        nElem.append(mesh['nElem']);
        pre_calculation(mesh);
        read_state(mesh,'..\\state_after_iteration\\alpha_%d_ada_iter_%d.txt' %(alpha,i));
        cl = find_coefs(mesh,U_inf)['cl'];
        cd = find_coefs(mesh,U_inf)['cd'];
        cm = find_coefs(mesh,U_inf)['cm'];
        coefs.append([cl,cd,cm]);
        p_versus_theta = plot_pressure_coef(mesh,folder,'adaptated_mesh_iter_%d_cp_versus_theta' %i,U_inf);
        coef_p.append(p_versus_theta);
        plot_mach(mesh,folder,'adaptated_mesh_iter_%d_Mach_number' %i);
    
    f1 = plt.figure(figsize=(8,6));
    l1,l2,l3 = plt.plot(nElem,coefs,'x-');
    plt.xlabel('Number of cells',fontsize =16);
    plt.ylabel('Coefficients',fontsize = 16);
    plt.grid();
    plt.legend([l1,l2,l3],['$c_l$','$c_d$','$c_m$'],prop={'size': 16});
    plt.savefig('..\\%s\\coefs_vs_num_cells.pdf' %folder,dpi=150);
    #plt.close(f1);
    f2 = plt.figure(figsize=(8,6));
    plt.plot(coef_p[1][:,0],coef_p[1][:,1],'x-',label = '$c_p$ after the $1^{st}$ adapative iteration');    
    plt.plot(coef_p[3][:,0],coef_p[3][:,1],'x-',label = '$c_p$ after the $3^{rd}$ adapative iteration');    
    plt.plot(coef_p[5][:,0],coef_p[5][:,1],'x-',label = '$c_p$ after the $5^{th}$ adapative iteration');    
    plt.xlabel(r"${\theta}$ deg",fontsize =16);
    plt.ylabel('Pressure coefficient $c_p$',fontsize =16);
    plt.grid();
    plt.legend();
    plt.savefig('..\\%s\\cp_vs_theta.pdf' %folder,dpi=150);
    plt.close(f2);
    return 0;

if __name__=="__main__":
    main()