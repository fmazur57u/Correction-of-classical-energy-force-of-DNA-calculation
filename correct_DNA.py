#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 10:51:41 2020

@author: Florian
"""

import numpy as np
from MDPlus.core import *
from scipy.interpolate import interp2d

##pdb file 182D with the intercalating dye molecule removed.
pdb = '182D_noInterc.pdb'

def correct_dna(pdb):

    a = open(pdb, 'r')
    f_list = []
    information_list_crd = []
    information_list_pdb = []
    crd_list = []
    pdb_list = []     
    del_list_start  = [0, 1, 2, 3, 4, 45, 46, 47, 48, 49, 50, 51, 52, 93, 94]
    del_list_middle = [0, 1, 2, 3, 4, 5,   6, 47, 48, 49, 50, 51, 52, 53, 54, 95]
    del_list_last   = [0, 1, 2, 3, 4, 5,   6, 47, 48, 49, 50, 51, 52, 53, 94]
    cr_start_list   = []
    cr_middle_list  = []
    cr_last_list    = []
    crd_start       = []
    pdb_start       = []
    crd_last        = []
    pdb_last        = []

    #Cut dna for obtain dyads and delete some atom of base

    ##only use generic base atoms: limited training set.
    dont_use_atoms = {'N9', 'C8', 'N7', 'H8', 'O6', \
                     'H1', 'N2', 'H21', 'H22', 'O2', \
                     'H6', 'H5', 'N4',  'H41', 'H42',\
                     'H2', 'N6', 'H61', 'H62', 'H3', \
                     'O4', 'C7', 'H71', 'H72', 'H73'}

    for line in a:
        as_Liste = line.split(" ")
        as_List = [elem for elem in as_Liste if elem.strip()]
        if as_List[0] == "ATOM":
            if as_List[2] not in dont_use_atoms: 

               
                ##wrong assumption here that the pdb lines are space separated, see the spec at
                ## pdb.org
                ## :  https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
                information_list_crd.append([float(as_List[len(as_List)-7]), as_List[len(as_List)-9],\
                                     as_List[len(as_List)-10], float(as_List[len(as_List)-6]),\
                               float(as_List[len(as_List)-5]), float(as_List[len(as_List)-4])])

                information_list_pdb.append(line)


    last_information = information_list_crd[-1]
    max_index = last_information[0]


    ##not clear what is happening here
    #Loop for create a dyads_list for crd format and for pdb format
    for p in range(int(2), int((max_index/2) + 2), int(2)):
        residue_list_crd = []
        residue_list_pdb = []
        for i in information_list_crd:
            if i[0] == p - 1 or\
               i[0] == p     or\
               i[0] == (max_index + 1) - p or\
               i[0] == (max_index + 1) - (p - 1): 

                residue_list_crd.append([i[3], i[4], i[5]])

        crd_list.append(residue_list_crd)
        for i in information_list_pdb:
            as_Liste = i.split(" ")
            as_List = [elem for elem in as_Liste if elem.strip()]
            if float(as_List[5]) == p -1 or \
               float(as_List[5]) == p or \
               float(as_List[5]) == (max_index + 1) - p or \
               float(as_List[5]) == (max_index + 1) - (p - 1):
                    residue_list_pdb.append(i)
        pdb_list.append(residue_list_pdb)

    #Delete some other atom in backbone for obtain a dyads with same number of atom
    for cr_start in range(len(crd_list[0])):
        cr_start_list.append(cr_start)

    ## If I remove the zero'th item in a list and then the first item,
    ## is that the same as removing the first item and then the zeroth item?
    ##
    ## It looks as if you have made this mistake here.
    for i_start_del in del_list_start:
        cr_start_list.remove(i_start_del)

    for i_start in cr_start_list:
        crd_start.append(crd_list[0][i_start])
        pdb_start.append(pdb_list[0][i_start])   

    np.savetxt('dyads_dna[0].pdb', pdb_start, fmt='%s')
    np.savetxt('dyads_dna[0].crd', crd_start)
    f0 = Fasu('dyads_dna[0].pdb', 'dyads_dna[0].crd')
    f_list.append(f0)
    for cr_middle in range(len(crd_list[1])):
        cr_middle_list.append(cr_middle)
    
    for i_middle_del in del_list_middle:
        cr_middle_list.remove(i_middle_del)
    
    for i_dyads_middle in range(1, len(crd_list) - 1):
        crd_middle = []
        pdb_middle = []
        for i_middle in cr_middle_list:
            crd_middle.append(crd_list[i_dyads_middle][i_middle])
            pdb_middle.append(pdb_list[i_dyads_middle][i_middle])
        np.savetxt('dyads_dna_%f.pdb' % i_dyads_middle, pdb_middle, fmt='%s')
        np.savetxt('dyads_dna_%f.crd' % i_dyads_middle, crd_middle)
        f_m = Fasu('dyads_dna_%f.pdb' % i_dyads_middle, 'dyads_dna_%f.crd' % i_dyads_middle)
        f_list.append(f_m)
    
    for cr_last in range(len(crd_list[-1])):
        cr_last_list.append(cr_last)   

    for i_last_del in del_list_last:
        cr_last_list.remove(i_last_del)

    for i_last in cr_last_list:
        crd_last.append(crd_list[-1][i_last])
        pdb_last.append(pdb_list[-1][i_last])
    np.savetxt('dyads_dna[-1].pdb', pdb_last, fmt='%s')
    np.savetxt('dyads_dna[-1].crd', crd_last)
    f_last = Fasu('dyads_dna[-1].pdb', 'dyads_dna[-1].crd')
    f_list.append(f_last)
    e0_space = np.loadtxt('e0_space_common.txt')
    e1_space = np.loadtxt('e1_space_common.txt')
    correct_E = []
    correct_F = []
    #Calulation of Fasu of reference
    ref_fasu = Fasu('dyads_GC_0.000000.pdb', 'traj_0.000000.crd')
    #loading of coordinates matrix and matrix of smooth surface
    coords = np.loadtxt('PCA_list_common.txt')
    E_matrix1 = np.loadtxt('E_matrix_e0_e1_common.txt')
    #Creation of interpolation spline for energy and force    
    spline = interp2d(e0_space, e1_space, E_matrix1, kind = 'cubic')
    derivative1 = interp2d.__call__(spline, e0_space, e1_space, dx=1, dy=0)
    derivative2 = interp2d.__call__(spline, e0_space, e1_space, dx=0, dy=1)
    gradient = interp2d(e0_space, e1_space, derivative1 + derivative2, kind = 'cubic') 
    #calculation of Cofasu of study geometry
    for i_fasu in range(len(f_list)):
        coords_PCA = []
        traj_data = []
        traj = Cofasu([ref_fasu, f_list[i_fasu]])
        traj.align()
        
        #addition of study geometry in coordinate matrix
        for i in range(len(traj[1])):
            for j in range(len(traj[1, i])):
                traj_data.append(traj[1, i, j])
        
        coords_PCA.append(traj_data)
        #Creation of coordinate matrix
        for i in range(len(coords)):
            coords_PCA.append(coords[i])
        coords_PCA = np.array(coords_PCA)
        #PCA
        print("pca, shape:", np.shape(coords_PCA))
        mean_state = np.sum(coords_PCA, axis=0)/len(coords_PCA[:,0])
        print("mean state, shape:", np.shape(mean_state))

        ##centre each configuration
        coords_PCA -= mean_state

        ##covariance matrix, no need to norm
        C = np.matmul(coords_PCA.transpose(), coords_PCA)

        print("covariance matrix, shape: ", np.shape(C))

        ##diagonalise Hermitian matrix
        eigvals, eigvecs = np.linalg.eigh(C)

        ##sort by eigvals, descending
        order   = np.argsort( eigvals )[::-1]
        eigvecs = eigvecs[:,order]
        eigvals = eigvals[order]
        #Valor of coordinates in two first principal component of geometry 
        e0 = np.dot(coords_PCA[0,:], eigvecs[:,0])
        e1 = np.dot(coords_PCA[0,:], eigvecs[:,1])           

 
        #Calculation of force and energy correction
        force_correct = gradient(e0, e1)[0]
        correct_energy = spline(e0, e1)[0]
        correct_E.append(correct_energy)
        correct_F.append(force_correct)
    #The dna is divided by 3 part
    part_E_1 = []
    part_E_2 = []
    part_E_3 = []
    part_F_1 = []
    part_F_2 = []
    part_F_3 = []
    for E in range(int(len(correct_E) / 3)):
        part_E_1.append(correct_E[E])
        part_F_1.append(correct_F[E])
    E_part1 = np.sum(part_E_1)
    F_part1 = np.sum(part_F_1)
    for EE in range(int(len(correct_E) / 3), int(2*len(correct_E) / 3)):
        part_E_2.append(correct_E[EE])
        part_F_2.append(correct_F[EE])
    E_part2 = np.sum(part_E_2)
    F_part2 = np.sum(part_F_2)
    for EEE in range(int(2*len(correct_E) / 3), int(len(correct_E))):
        part_E_3.append(correct_E[EEE])
        part_F_3.append(correct_F[EEE])
    E_part3 = np.sum(part_E_3)
    F_part3 = np.sum(part_F_3)
    #Correction with division by 2 with the all part
    Correction_E_dna = E_part1/2 + E_part2/2 + E_part3/2
    Correction_F_dna = F_part1/2 + F_part2/2 + F_part3/2
    #Correction[0] = energy correction and Correction[1] = force correction
    Correction = [Correction_E_dna, Correction_F_dna]
    return Correction
     
##test code, does not run when this code is imported as a module, only when 
##called as an appli   
if __name__ == "__main__":  
     print(correct_dna(pdb))

