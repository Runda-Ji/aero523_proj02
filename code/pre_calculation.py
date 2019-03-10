# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 10:51:26 2018

@author: rundaji
"""

from my_class import EDGE
from scipy.sparse import lil_matrix
import numpy as np

def pre_calculation(mesh):
    conn_vertex = find_adjacent_cell(mesh);
    furnish_boundary(mesh,conn_vertex);
    furnish_edge(mesh);
    return 0;

#--------------------------------------------------------------------------

def find_adjacent_cell(mesh):
    nNode = mesh['nNode'];
    nElem = mesh['nElem'];
    Elem = mesh['Elems'];
    conn_vertex = lil_matrix((nNode,nNode),dtype = int);
    # we want to store the global edge # in sparse matrix conn_vertex
    # however the index has to start from 1
    for i in range(0,nElem):
        v_0,v_1,v_2 = Elem[i].vertex;
        e_0,e_1,e_2 = Elem[i].edge;
        check(e_0,conn_vertex,mesh); #check if e_0 already in the list
        check(e_1,conn_vertex,mesh); #check if e_1 already in the list
        check(e_2,conn_vertex,mesh); #check if e_2 already in the list
    return conn_vertex;

def check(e_curr,conn_vertex,mesh):
    #check if the edge is already in the v^{th} line
    n0 = e_curr['n0'];
    n1 = e_curr['n1'];
    t = e_curr['t'];
    e = e_curr['e'];
    if conn_vertex[n0,n1] == 0:
        #if edge does NOT exist in matrix, add this edge to the matrix
        t1 = t; e1 = e;
        global_edge_no = mesh['nEdge'];
        mesh['Edges'].append(EDGE(n0,n1,t1,e1,'No exterior cell','Boundary condition',None,None,None));
        mesh['nEdge'] = mesh['nEdge'] + 1;
        conn_vertex[n0,n1] = global_edge_no + 1;
        conn_vertex[n1,n0] = global_edge_no + 1;
        mesh['Elems'][t].edge[e] = global_edge_no;
        #the index stored in the sparse matrix is index of edge + 1
    else:
        #if the edge already exist in the list
        t2 = t; e2 = e;
        #decide which edge alread in the list
        global_edge_no = conn_vertex[n0,n1]-1;
        mesh['Elems'][t2].edge[e2] = global_edge_no;
        t1 = mesh['Edges'][global_edge_no].t1;
        e1 = mesh['Edges'][global_edge_no].e1;
        #store the info of new edge
        mesh['Edges'][global_edge_no].t2 = t2;
        mesh['Edges'][global_edge_no].e2 = e2;
        create_adjacent_cell(t1,e1,t2,e2,mesh);
    return 0;
    
def create_adjacent_cell(t1,e1,t2,e2,mesh):
    adj_info_1 = {'loc_edge': e1, 'adj_edge': e2, 'adj_tri': t2};
    mesh['Elems'][t1].adj_cell.append(adj_info_1);
    adj_info_2 = {'loc_edge': e2, 'adj_edge': e1, 'adj_tri': t1};
    mesh['Elems'][t2].adj_cell.append(adj_info_2);
    return 0;

#--------------------------------------------------------------------------

def furnish_boundary(mesh,conn_vertex):
    for i in range(0,mesh['nBGroup']):
        #5 boundary groups
        for j in range(0,mesh['boundary'][i].nBFace):
            Title = mesh['boundary'][i].Title;
            n0 = mesh['boundary'][i].B[j]['n0'];
            n1 = mesh['boundary'][i].B[j]['n1'];
            global_edge_no = conn_vertex[n0,n1] - 1;
            mesh['Edges'][global_edge_no].e2 = Title;
            mesh['boundary'][i].B[j] = global_edge_no;
    return 0;

#--------------------------------------------------------------------------

def furnish_edge(mesh):
    for i in range(0,mesh['nEdge']):
        A = mesh['Edges'][i].node_0;
        B = mesh['Edges'][i].node_1;
        
        xA,yA = mesh['node_pos'][A];
        xB,yB = mesh['node_pos'][B];
        
        #find length
        vec_len = np.sqrt((xA - xB)**2 + (yA - yB)**2);
        mesh['Edges'][i].length = vec_len;
        
        #find normal vector
        vec = np.array([(yB-yA), (xA-xB)]);
        mesh['Edges'][i].norm_vec = vec/vec_len;
    return 0;

