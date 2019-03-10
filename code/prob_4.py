# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 10:23:21 2018

@author: rundaji
"""

#from my_class import EDGE, CELL, BOUNDARY

from read_data import read_data
from pre_calculation import pre_calculation
from initialization import initialization,find_U_inf,read_state
from time_iteration import time_iteration
from adaptation import adaptation
from post_processing import plot_mesh,post_processing,plot_color_flag

def main():
    global mesh;
    alpha = 5;
    folder = ('figure_prob_4\\alpha_%d' %alpha);
    mesh = read_data('..\\mesh\\capsule.gri');
    plot_mesh(mesh,folder,'Base_line_mesh');
    pre_calculation(mesh);
    #main iteration
    initialization(mesh,alpha);
    U_inf = find_U_inf(alpha);
    n_Time_Step = 5000;
    history = time_iteration(mesh,n_Time_Step,U_inf,'Base_line');
    #post processing
    Mach = post_processing(mesh,history,folder,'Base_line',U_inf);
    
    for i in range(1,6):
        # create the ith adaptation
        color_flag = adaptation(mesh,Mach,alpha,i);
        plot_color_flag(mesh,color_flag,folder,'adaptated_mesh_iter_%d' %i)
        # read in the mesh after ith adaptation
        mesh = read_data('..\\mesh\\alpha_%d_ada_iter_%d.gri' %(alpha,i));
        plot_mesh(mesh,folder,'adaptated_mesh_iter_%d_mesh' %i);
        read_state(mesh,'..\\state_before_iteration\\alpha_%d_ada_iter_%d.txt' %(alpha,i));
        pre_calculation(mesh);
        #run simulation on the ith adaptative mesh
        history = time_iteration(mesh,n_Time_Step,U_inf, 'alpha_%d_ada_iter_%d' %(alpha,i));
        Mach = post_processing(mesh,history,folder,'adaptated_mesh_iter_%d' %i,U_inf);
    return 0;

if __name__=="__main__":
    main()