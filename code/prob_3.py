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
from post_processing import post_processing
#------------------------------------------------------------------------------
from post_processing import plot_mach

def main():
    global mesh;
    alpha = 5;
    folder = ('figure_prob_3\\');
    mesh = read_data('..\\mesh\\capsule.gri');
    pre_calculation(mesh);
    #main iteration
    initialization(mesh,alpha);
    U_inf = find_U_inf(alpha);
    n_Time_Step = 5000;
    history = time_iteration(mesh,n_Time_Step,U_inf,'Base_line');
    #post processing
    Mach = post_processing(mesh,history,folder,'Base_line',U_inf);
    return 0;

if __name__=="__main__":
    main()