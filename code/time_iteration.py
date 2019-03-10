# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 23:28:43 2018

@author: Runda Ji
"""

import numpy as np
from flux import flux, HLLE_flux, wall_flux

def time_iteration(mesh,nTime_step,U_inf,fname):
    L2_error = [];
    cl = [];
    cd = [];
    cm = [];
    for i in range(0,nTime_step):
        refresh_R(mesh,U_inf);
        compute_time_step(mesh,0.5);
        Res = single_time_step(mesh);
        err = np.sqrt(np.sum(Res));
        L2_error.append(err);
        coefs = find_coefs(mesh,U_inf);
        cl.append(coefs['cl']);
        cd.append(coefs['cd']);
        cm.append(coefs['cm']);
        if L2_error[i] < 1e-5:
            auto_save(mesh,'..\\state_after_iteration\\%s.txt' %fname);
            save_history(L2_error,cl,cd,cm,'..\\history_after_iteration\\%s.txt' %fname);
            break;
        if np.remainder(i,100) == 0:
            print('iteration=%d,L2_error=%f' %(i,L2_error[i]));
            auto_save(mesh,'..\\state_auto_save\\%s_time_iter_%d.txt' %(fname,i));
    history = {'L2_error':L2_error,'cl':cl, 'cd':cd, 'cm':cm};
    return history;

#--------------------------------------------------------------------------

def refresh_R(mesh,U_inf):
    #empty R for all cells
    for i in range(0,mesh['nElem']):
        mesh['Elems'][i].R = np.array([0.0,0.0,0.0,0.0]);
    for i in range(0,mesh['nEdge']):
        edge = mesh['Edges'][i];
        n = edge.norm_vec;
        delta_l = edge.length;
        if type(edge.e2) == type(1):
            #for interior cells, e2 is an integer
            t1 = edge.t1;
            t2 = edge.t2;
            U_L = mesh['Elems'][t1].state;
            U_R = mesh['Elems'][t2].state;
            F,s = HLLE_flux(U_L,U_R,n);
            
            #-----DEBUG-----DEBUG-----DEBUG-----DEBUG-----DEBUG-----DEBUG------
            if s == -1:
                print('Negative pressure @ the %d^th edge' %i);
                n0 = edge.node_0;
                n1 = edge.node_1;
                n0_pos = mesh['node_pos'][n0];
                n1_pos = mesh['node_pos'][n1];
                print('node_0 = %d, node_1 = %d' %(n0,n1));
                print('node_0 position %f %f' %(n0_pos[0],n0_pos[1]));
                print('node_1 position %f %f' %(n1_pos[0],n1_pos[1]));
            #-----DEBUG-----DEBUG-----DEBUG-----DEBUG-----DEBUG-----DEBUG------
            
            mesh['Elems'][t1].R = mesh['Elems'][t1].R + F*delta_l;
            mesh['Elems'][t2].R = mesh['Elems'][t2].R - F*delta_l;
            mesh['Edges'][i].s = s;  
        if edge.e2 == 'Bottom' or edge.e2 == 'Top' or edge.e2 == 'Left':
            #for full state
            t = edge.t1;
            U_L = mesh['Elems'][t].state;
            U_R = U_inf;
            F,s = HLLE_flux(U_L,U_R,n);
            mesh['Elems'][t].R = mesh['Elems'][t].R + F*delta_l;
            mesh['Edges'][i].s = s;    
        if edge.e2 == 'Capsule':
            #for wall
            t = edge.t1;
            U_L = mesh['Elems'][t].state;
            F,s = wall_flux(U_L,n);
            
            #-----DEBUG-----DEBUG-----DEBUG-----DEBUG-----DEBUG-----DEBUG------
            if s == -1:
                print('Negative pressure @ the %d^th edge' %i);
                n0 = edge.node_0;
                n1 = edge.node_1;
                n0_pos = mesh['node_pos'][n0];
                n1_pos = mesh['node_pos'][n1];
                print('node_0 = %d, node_1 = %d' %(n0,n1));
                print('node_0 position %f %f' %(n0_pos[0],n0_pos[1]));
                print('node_1 position %f %f' %(n1_pos[0],n1_pos[1]));
            #-----DEBUG-----DEBUG-----DEBUG-----DEBUG-----DEBUG-----DEBUG------
            
            mesh['Elems'][t].R = mesh['Elems'][t].R + F*delta_l;
            mesh['Edges'][i].s = s;   
        if edge.e2 == 'Right':
            #for downstream
            t = edge.t1;
            U_L = mesh['Elems'][t].state;
            prop = flux(U_L);
            F = prop['F'][0]*n[0] + prop['F'][1]*n[1];
            c = prop['a'];
            Vn = np.dot(prop['V'],n);
            mesh['Elems'][t].R = mesh['Elems'][t].R + F*delta_l;
            mesh['Edges'][i].s = np.abs(Vn) + c;      
    return 0;

#--------------------------------------------------------------------------

def compute_time_step(mesh,CFL):
    for i in range(0,mesh['nElem']):
        SUM = 0;
        for j in range(0,3):
            global_edge_no = mesh['Elems'][i].edge[j];
            s = mesh['Edges'][global_edge_no].s;
            delta_l = mesh['Edges'][global_edge_no].length;
            SUM = SUM + s*delta_l;
        mesh['Elems'][i].dt = 2*CFL/SUM;
    return 0;
    
def single_time_step(mesh):
    Residual = 0;
    for i in range(0,mesh['nElem']):
        dt = mesh['Elems'][i].dt;
        R = mesh['Elems'][i].R;
        Residual = Residual + R**2;
        mesh['Elems'][i].state = mesh['Elems'][i].state - dt*R;
    return Residual;

#--------------------------------------------------------------------------

def auto_save(mesh,fname):
    file = open(fname, 'w');
    for i in range(0,mesh['nElem']):
        state = mesh['Elems'][i].state;
        file.write('%f %f %f %f\n' %(state[0],state[1],state[2],state[3]));
    file.close();
    return 0;

#--------------------------------------------------------------------------
def save_history(L2_error,cl,cd,cm,fname):
    file = open(fname, 'w');
    for i in range(0,len(L2_error)):
        file.write('%f %f %f %f\n' %(L2_error[i],cl[i],cd[i],cm[i]));
    file.close();
    return 0;    


#--------------------------------------------------------------------------

def find_coefs(mesh,U_inf):
    gamma = 1.3;
    rho_inf = U_inf[0];
    u_inf = U_inf[1]/U_inf[0];
    v_inf = U_inf[2]/U_inf[0];
    d = 0.6*2;
    norm_factor_force = 0.5*rho_inf*(u_inf**2 + v_inf**2)*d;
    norm_factor_momentum = 0.5*rho_inf*(u_inf**2 + v_inf**2)*d**2;
    for i in range(0,mesh['nBGroup']):
        if mesh['boundary'][i].Title == 'Capsule':
            capsule = i;
            break;
    F = np.array([0,0]);
    M = 0;
    for i in range(0,mesh['boundary'][capsule].nBFace):
        global_edge_no = mesh['boundary'][capsule].B[i];
        edge = mesh['Edges'][global_edge_no];
        state = mesh['Elems'][edge.t1].state;
        rho = state[0];
        u = state[1]/state[0];
        v = state[2]/state[0];
        E = state[3]/state[0];
        p = (gamma-1)*(rho*E-0.5*rho*(u**2+v**2));
        F = F + p*edge.norm_vec*edge.length;
        pos_n0 = mesh['node_pos'][edge.node_0];
        pos_n1 = mesh['node_pos'][edge.node_1];
        r = 0.5*(pos_n0 + pos_n1);
        M = M + p*np.cross(r,edge.norm_vec)*edge.length;
    D = F[0]; L = F[1];
    cl = L/norm_factor_force;
    cd = D/norm_factor_force;
    cm = M/norm_factor_momentum;
    coefs = {'cl':cl,'cd':cd,'cm':cm};
    return coefs;