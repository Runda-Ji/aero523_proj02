# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 19:34:36 2018

@author: Runda Ji
"""

import numpy as np;
import matplotlib.pyplot as plt
import operator

def post_processing(mesh,history,folder,name, U_inf):
    plot_history(history,folder,'%s_history' %name);
    Mach = plot_mach(mesh,folder,'%s_Mach_number' %name);
    plot_pressure_coef(mesh,folder,'%s_cp_versus_theta' %name, U_inf)
    return Mach;

def plot_mesh(mesh,folder,title):
    x = mesh['node_pos'][:,0];
    y = mesh['node_pos'][:,1];
    tri = [None]*mesh['nElem'];
    for i in range(0,mesh['nElem']):
        tri[i] = mesh['Elems'][i].vertex;
    tri = np.asarray(tri);
    f1 = plt.figure(figsize=(8,6));
    plt.triplot(x, y, tri, 'k-', lw=0.5);
    plt.axis('equal');
    plt.savefig('..\\%s\\%s.pdf' %(folder,title),dpi=150);
    plt.close(f1);
    f2 = plt.figure(figsize=(8,6));
    plt.triplot(x, y, tri, 'k-', lw=0.5);
    plt.axis('equal');
    plt.xlim((-3, 2));
    plt.ylim((-3, 3));
    plt.savefig('..\\%s\\%s_zoom_in.pdf' %(folder,title),dpi=150);
    plt.close(f2);
    return 0;
    
def plot_mach(mesh,folder,title):
    gamma = 1.3;
    x = mesh['node_pos'][:,0];
    y = mesh['node_pos'][:,1];
    tri = [None]*mesh['nElem'];
    M = [None]*mesh['nElem'];
    for i in range(0,mesh['nElem']):
        tri[i] = mesh['Elems'][i].vertex;
        state = mesh['Elems'][i].state;
        rho = state[0];
        u = state[1]/state[0];
        v = state[2]/state[0];
        E = state[3]/state[0];
        p = (gamma-1)*(rho*E-0.5*rho*(u**2+v**2));
        a = np.sqrt(gamma*p/rho);
        M[i] = np.sqrt(u**2+v**2)/a;
    tri = np.asarray(tri);
    f1 = plt.figure(figsize=(10,6));
    plt.tripcolor(x, y, tri, M, edgecolors='k', cmap=plt.cm.jet);
    plt.axis('equal');
    plt.colorbar();
    plt.savefig('..\\%s\\%s.pdf' %(folder,title),dpi=150);
    plt.close(f1);
    #--------------------------------------------------------------------------
    f2 = plt.figure(figsize=(10,6));
    plt.tripcolor(x, y, tri, M, edgecolors='k', cmap=plt.cm.jet);
    plt.axis('equal');
    plt.colorbar();
    plt.xlim((-3, 4));
    plt.ylim((-3, 3));
    plt.savefig('..\\%s\\%s_zoom_in.pdf' %(folder,title),dpi=150);
    plt.close(f2);
    return M;

def plot_history(history,folder,fname):
    L2 = history['L2_error'];
    f1 = plt.figure(figsize=(8,6));
    plt.semilogy(L2,'k-');
    plt.xlabel('Iteration',fontsize =16);
    plt.ylabel('$L_2$ error',fontsize =16);
    plt.grid();
    plt.savefig('..\\%s\\%s_L2_error.pdf' %(folder,fname),dpi=150);
    plt.close(f1);
    #--------------------------------------------------------------------------
    f2 = plt.figure(figsize=(8,6));
    cl = history['cl'];
    cd = history['cd'];
    cm = history['cm'];
    l1 = plt.plot(cl,label='$c_l$',lw=1.5);
    l2 = plt.plot(cd,label='$c_d$',lw=1.5);
    l3 = plt.plot(cm,label='$c_m$',lw=1.5);
    plt.xlabel('Iteration',fontsize =16);
    plt.ylabel('Coefficients',fontsize = 16);
    plt.grid();
    plt.legend(prop={'size': 16});
    plt.savefig('..\\%s\\%s_coefficients.pdf' %(folder,fname),dpi=150);
    plt.close(f2);
    return 0;

def plot_color_flag(mesh,color_flag,folder,title):
    x = mesh['node_pos'][:,0];
    y = mesh['node_pos'][:,1];
    tri = [None]*mesh['nElem'];
    for i in range(0,mesh['nElem']):
        tri[i] = mesh['Elems'][i].vertex;
    tri = np.asarray(tri);
    f1 = plt.figure(figsize=(10,6));
    plt.tripcolor(x, y, tri, color_flag, edgecolors='k', cmap=plt.cm.jet);
    plt.axis('equal');
    plt.colorbar();
    plt.savefig('..\\%s\\%s_color_flag.pdf' %(folder,title),dpi=150);
    plt.close(f1);
    #--------------------------------------------------------------------------
    f2 = plt.figure(figsize=(10,6));
    plt.tripcolor(x, y, tri, color_flag, edgecolors='k', cmap=plt.cm.jet);
    plt.axis('equal');
    plt.colorbar();
    plt.xlim((-3, 2));
    plt.ylim((-3, 3));
    plt.savefig('..\\%s\\%s_color_flag_zoom_in.pdf' %(folder,title),dpi=150);
    plt.close(f2);
    return 0;

def plot_pressure_coef(mesh,folder,fname,U_inf):
    gamma = 1.3;
    rho_inf = U_inf[0];
    u_inf = U_inf[1]/U_inf[0];
    v_inf = U_inf[2]/U_inf[0];
    E_inf =  U_inf[3]/U_inf[0];
    p_inf = (gamma-1)*(rho_inf*E_inf-0.5*rho_inf*(u_inf**2+v_inf**2));
    n_x = np.array([1,0]);
    c = np.array([0.8,0]);
    norm_factor_pressure = 0.5*rho_inf*(u_inf**2 + v_inf**2);
    for i in range(0,mesh['nBGroup']):
        if mesh['boundary'][i].Title == 'Capsule':
            capsule = i;
            break;
    p_versus_theta = [];
    for i in range(0,mesh['boundary'][capsule].nBFace):
        global_edge_no = mesh['boundary'][capsule].B[i];
        edge = mesh['Edges'][global_edge_no];
        pos_n0 = mesh['node_pos'][edge.node_0];
        pos_n1 = mesh['node_pos'][edge.node_1];
        if pos_n0[0] <= 0 and pos_n1[0] <= 0:
            mid_pt = 0.5*(pos_n0 + pos_n1);
            state = mesh['Elems'][edge.t1].state;
            rho = state[0];
            u = state[1]/state[0];
            v = state[2]/state[0];
            E = state[3]/state[0];
            p = (gamma-1)*(rho*E-0.5*rho*(u**2+v**2));
            cp = (p-p_inf)/norm_factor_pressure;
            vec = mid_pt - c;
            len_vec = np.sqrt(np.sum(vec**2));
            cos_theta = np.dot(vec,-n_x)/len_vec;
            theta = np.arccos(cos_theta)*180/np.pi;
            if mid_pt[1] < 0:
                theta = -theta;
            p_versus_theta.append([theta,cp]);
    p_versus_theta.sort(key = operator.itemgetter(0),reverse = False);
    #--------------------------------------------------------------------------
    f = plt.figure(figsize=(8,6));
    p_versus_theta = np.asarray(p_versus_theta);
    x = p_versus_theta[:,0];
    y = p_versus_theta[:,1];
    plt.plot(x,y,'kx-');
    plt.xlabel(r"${\theta}$ deg",fontsize =16);
    plt.ylabel('Pressure coefficient $c_p$',fontsize =16);
    plt.grid();
    plt.savefig('..\\%s\\%s.pdf' %(folder,fname),dpi=150);
    plt.close(f);
    return p_versus_theta;