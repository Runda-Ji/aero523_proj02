# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 10:45:58 2018

@author: rundaji
"""

from my_class import EDGE, CELL
import numpy as np
import operator
from time_iteration import auto_save

def adaptation(mesh,Mach,alpha,ada_iter):
    flag_edge = error_indicator(mesh,Mach);
    new_nodes,new_edges = split_edge(mesh,flag_edge);
    color_flag = split_elem(mesh,flag_edge,new_nodes);
    split_boundary(mesh,flag_edge,new_edges);
    write_gri_file(mesh,alpha,ada_iter);
    auto_save(mesh,'..\\state_before_iteration\\alpha_%d_ada_iter_%d.txt' %(alpha,ada_iter));
    return color_flag;

def error_indicator(mesh,Mach):
    #find the eps_e
    gamma = 1.3;
    eps_e = [None]*mesh['nEdge'];
    for i in range(0,mesh['nEdge']):
        edge = mesh['Edges'][i];
        if type(edge.e2) == type(1):
            #for interior cells, e2 is an integer
            Mkplus = Mach[edge.t1];
            Mkminus = Mach[edge.t2];
            eps_e[i] = [i,np.abs(Mkplus-Mkminus)*edge.length];
        if edge.e2 == 'Bottom' or edge.e2 == 'Top' or edge.e2 == 'Left' or edge.e2 == 'Right':
            eps_e[i] = [i,0];
        if edge.e2 == 'Capsule':
            n = edge.norm_vec;
            elem = mesh['Elems'][edge.t1];
            rho = elem.state[0];
            u = elem.state[1]/elem.state[0];
            v = elem.state[2]/elem.state[0];
            E = elem.state[3]/elem.state[0];
            p = (gamma-1)*(rho*E-0.5*rho*(u**2+v**2));
            a = np.sqrt(gamma*p/rho);
            V = np.array([u,v]);
            Vn = np.dot(V,n);
            Mkvertical = np.sqrt(np.sum(Vn**2))/a;
            eps_e[i] = [i,np.abs(Mkvertical)*edge.length];    
    #sort and flag the edges
    f = 0.03;
    eps_e.sort(key = operator.itemgetter(1),reverse = True);
    flag_edge = np.zeros(mesh['nEdge']);
    #step (1) flag a small fraction of edges
    for i in range(0,round(f*mesh['nEdge'])):
        global_edge_no = eps_e[i][0];
        flag_edge[global_edge_no] = 1;
    #step (2) smooth the flag
    for i in range(0,mesh['nElem']):
        f = 0;
        #check if any edge of cell i is flagged
        for j in range(0,3):
            global_edge_no = mesh['Elems'][i].edge[j];
            if flag_edge[global_edge_no] == 1:
                f = 1;
        #if so then flag all the edges
        if f == 1:
            for j in range(0,3):
                global_edge_no = mesh['Elems'][i].edge[j];
                flag_edge[global_edge_no] = 1;
    return flag_edge;

def split_edge(mesh,flag_edge):    
    center = np.array([0.8,0]);
    R = 1;
    nEdge_old = mesh['nEdge'];
    new_nodes = [None]*nEdge_old;
    new_edges = [None]*nEdge_old;
    for i in range(0,nEdge_old):
        if flag_edge[i] == 1:
            edge = mesh['Edges'][i];
            x0,y0 = mesh['node_pos'][edge.node_0];
            x1,y1 = mesh['node_pos'][edge.node_1];
            midpt_pos = np.array([0.5*(x0+x1),0.5*(y0+y1)]);
            #correct the midpt position for heat shield
            if edge.e2 == 'Capsule':
                if x0 <= 0 and x1 <= 0:
                    vec = midpt_pos - center;
                    length = np.sqrt(np.sum(vec**2));
                    vec_corrected = (R/length)*vec;
                    midpt_pos = center + vec_corrected;
            #add the midpt to current list
            mesh['node_pos'] = np.append(mesh['node_pos'], [midpt_pos], axis = 0);
            new_node_no = mesh['nNode'];
            mesh['nNode'] = mesh['nNode'] + 1;
            new_nodes[i] = new_node_no;
            #split the edge            
            new_edge_0 = EDGE(edge.node_0, new_node_no, None, None, None, None, None, None, None);
            new_edge_1 = EDGE(new_node_no, edge.node_1, None, None, None, None, None, None, None);
            mesh['Edges'][i] = new_edge_0; # replace the current edge with edge0
            mesh['Edges'].append(new_edge_1); # append the edge1
            new_edge_no = mesh['nEdge'];
            mesh['nEdge'] = mesh['nEdge'] + 1;
            new_edges[i] = {'new_edge_0': i, 'new_edge_1': new_edge_no};
    return new_nodes,new_edges;

def split_elem(mesh,flag_edge,new_nodes):
    global edge_splitted;
    nElem_old = mesh['nElem'];
    color_flag = [0]*nElem_old;
    for i in range(0,nElem_old):
        #count the flagged edges of cell
        count = 0;
        edge_splitted = [];
        cell = mesh['Elems'][i];
        state = cell.state;
        for j in range(0,3):
            global_edge_no = cell.edge[j];
            if flag_edge[global_edge_no] == 1:
                count = count + 1;
                edge_splitted.append({'loc_edge': j, 'new_node_no': new_nodes[global_edge_no]});
        if count == 1:
            loc_v0 = edge_splitted[0]['loc_edge'];
            loc_v1 = np.remainder(loc_v0 + 1,3);
            loc_v2 = np.remainder(loc_v0 + 2,3);
            global_v0 = cell.vertex[loc_v0];
            global_v1 = cell.vertex[loc_v1];
            global_v2 = cell.vertex[loc_v2];
            new_node_no = edge_splitted[0]['new_node_no'];
            #adding the split line to list of edge
            new_edge = EDGE(global_v0, new_node_no, None, None, None, None, None, None, None);
            mesh['Edges'].append(new_edge);
            mesh['nEdge'] = mesh['nEdge'] + 1;
            #split the cell
            new_tri_0 = CELL([global_v0,global_v1,new_node_no], None, None, None, state, None, None);
            new_tri_1 = CELL([global_v0,new_node_no,global_v2], None, None, None, state, None, None)
            mesh['Elems'][i] = new_tri_0;
            color_flag[i] = 1;
            mesh['Elems'].append(new_tri_1);
            color_flag.append(1);
            mesh['nElem'] = mesh['nElem'] + 1;
        if count == 2:
            #always set the non-splitted edge as e0 corresponding to v0
            if edge_splitted[0]['loc_edge'] == 0 and edge_splitted[1]['loc_edge'] == 1:
                loc_v0 = 2; loc_v1 = 0; loc_v2 = 1;
                new_node_no_0 = edge_splitted[1]['new_node_no'];
                new_node_no_1 = edge_splitted[0]['new_node_no'];
            if edge_splitted[0]['loc_edge'] == 0 and edge_splitted[1]['loc_edge'] == 2:
                loc_v0 = 1; loc_v1 = 2; loc_v2 = 0;
                new_node_no_0 = edge_splitted[0]['new_node_no'];
                new_node_no_1 = edge_splitted[1]['new_node_no'];
            if edge_splitted[0]['loc_edge'] == 1 and edge_splitted[1]['loc_edge'] == 2:
                loc_v0 = 0; loc_v1 = 1; loc_v2 = 2;
                new_node_no_0 = edge_splitted[1]['new_node_no'];
                new_node_no_1 = edge_splitted[0]['new_node_no'];
            global_v0 = cell.vertex[loc_v0];
            global_v1 = cell.vertex[loc_v1];
            global_v2 = cell.vertex[loc_v2];
            #adding the 1^th split line to list of edge
            new_edge_0 = EDGE(new_node_no_0, new_node_no_1, None, None, None, None, None, None, None);
            mesh['Edges'].append(new_edge_0);
            mesh['nEdge'] = mesh['nEdge'] + 1;
            #adding the 1^th cell to the list of elements
            new_tri_0 = CELL([global_v0,new_node_no_0,new_node_no_1], None, None, None, state, None, None);
            mesh['Elems'][i] = new_tri_0;
            color_flag[i] = 2;
            #determing the 2^nd solit line
            larger = find_large_angle(mesh, global_v0,global_v1,global_v2);
            if larger == 1:
                new_edge_1 = EDGE(global_v1, new_node_no_1, None, None, None, None, None, None, None);
                new_tri_1 = CELL([new_node_no_0, global_v1, new_node_no_1], None, None, None, state, None, None);
                new_tri_2 = CELL([new_node_no_1, global_v1, global_v2], None, None, None, state, None, None);
            else:
                new_edge_1 = EDGE(global_v2, new_node_no_0, None, None, None, None, None, None, None);
                new_tri_1 = CELL([new_node_no_0, global_v1, global_v2], None, None, None, state, None, None);
                new_tri_2 = CELL([new_node_no_1, new_node_no_0, global_v2], None, None, None, state, None, None);
            mesh['Edges'].append(new_edge_1);
            mesh['nEdge'] = mesh['nEdge'] + 1;
            mesh['Elems'].append(new_tri_1);
            mesh['Elems'].append(new_tri_2);
            color_flag.append(2);
            color_flag.append(2);
            mesh['nElem'] = mesh['nElem'] + 2;
        if count == 3:
            loc_v0 = edge_splitted[0]['loc_edge']; #v0 = 0
            loc_v1 = edge_splitted[1]['loc_edge']; #v0 = 0
            loc_v2 = edge_splitted[2]['loc_edge']; #v0 = 0
            global_v0 = cell.vertex[loc_v0];
            global_v1 = cell.vertex[loc_v1];
            global_v2 = cell.vertex[loc_v2];
            new_node_no_0 = edge_splitted[0]['new_node_no'];
            new_node_no_1 = edge_splitted[1]['new_node_no'];
            new_node_no_2 = edge_splitted[2]['new_node_no'];
            new_edge_0 = EDGE(new_node_no_0, new_node_no_1, None, None, None, None, None, None, None);
            new_edge_1 = EDGE(new_node_no_1, new_node_no_2, None, None, None, None, None, None, None);
            new_edge_2 = EDGE(new_node_no_2, new_node_no_0, None, None, None, None, None, None, None);
            new_tri_0 = CELL([global_v0, new_node_no_2, new_node_no_1], None, None, None, state, None, None);
            new_tri_1 = CELL([new_node_no_2, global_v1, new_node_no_0], None, None, None, state, None, None);
            new_tri_2 = CELL([new_node_no_2, new_node_no_0, new_node_no_1], None, None, None, state, None, None);
            new_tri_3 = CELL([new_node_no_1, new_node_no_0, global_v2], None, None, None, state, None, None);
            mesh['Edges'].append(new_edge_0);
            mesh['Edges'].append(new_edge_1);
            mesh['Edges'].append(new_edge_2);
            mesh['nEdge'] = mesh['nEdge'] + 3;
            mesh['Elems'][i] = new_tri_0;
            color_flag[i] = 3;
            mesh['Elems'].append(new_tri_1);
            mesh['Elems'].append(new_tri_2);
            mesh['Elems'].append(new_tri_3);
            color_flag.append(3);
            color_flag.append(3);
            color_flag.append(3);
            mesh['nElem'] = mesh['nElem'] + 3;
    return color_flag;

def split_boundary(mesh,flag_edge,new_edges):
    for i in range(0,mesh['nBGroup']):
        nBFace_old = mesh['boundary'][i].nBFace;
        for j in range(0,nBFace_old):
            global_edge_no = mesh['boundary'][i].B[j];
            if flag_edge[global_edge_no] == 1:
                new_edge = new_edges[global_edge_no]['new_edge_1'];
                mesh['boundary'][i].B.append(new_edge);
                mesh['boundary'][i].nBFace = mesh['boundary'][i].nBFace + 1;
    return 0;
                
def find_large_angle(mesh, glbal_v0,glbal_v1,glbal_v2):
    #decide the large angle
    pos_v0 = mesh['node_pos'][glbal_v0];
    pos_v1 = mesh['node_pos'][glbal_v1];
    pos_v2 = mesh['node_pos'][glbal_v2];
    vec_20 = pos_v0 - pos_v2;
    len_20 = np.sqrt(np.sum(vec_20**2));
    vec_10 = pos_v0 - pos_v1;
    len_10 = np.sqrt(np.sum(vec_10**2));
    vec_12 = pos_v2 - pos_v1;
    len_12 = np.sqrt(np.sum(vec_12**2));
    cos_theta_1 = np.dot(vec_10,vec_12)/(len_10*len_12);
    theta_1 =  np.arccos(cos_theta_1);
    cos_theta_2 = np.dot(vec_20,-vec_12)/(len_20*len_12);
    theta_2 =  np.arccos(cos_theta_2);
    if theta_1 > theta_2:
        larger = 1;
    else:
        larger = 2;
    return larger;

def write_gri_file(mesh,alpha,ada_iter):
    file = open('..\\mesh\\alpha_%d_ada_iter_%d.gri' %(alpha,ada_iter), 'w');
    file.write('%d %d %d\n' %(mesh['nNode'],mesh['nElem'],2));
    for i in range(0,mesh['nNode']):
        pos = mesh['node_pos'][i];
        file.write('%f %f\n' %(pos[0], pos[1]));
    file.write('%d\n' %(mesh['nBGroup']));
    for i in range(0,mesh['nBGroup']):
        file.write('%d %d %s\n' %(mesh['boundary'][i].nBFace, mesh['boundary'][i].nf, mesh['boundary'][i].Title));
        for j in range(0,mesh['boundary'][i].nBFace):
            global_edge_no = mesh['boundary'][i].B[j];
            node_0 = mesh['Edges'][global_edge_no].node_0 + 1; #the index start from 1
            node_1 = mesh['Edges'][global_edge_no].node_1 + 1; #the index start from 1
            file.write('%d %d\n' %(node_0, node_1));
    file.write('%d %d %s\n' %(mesh['nElem'],1,'triangles'));
    for i in range(0,mesh['nElem']):
        elem = mesh['Elems'][i];
        v0,v1,v2 = elem.vertex;
        v0 = v0 + 1; v1 = v1 + 1; v2 = v2 + 1;
        file.write('%d %d %d\n' %(v0,v1,v2));
    file.close();
    return 0;