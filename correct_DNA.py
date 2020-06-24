# -*- coding: utf-8 -*-
"""
Created on Tue May 19 10:51:41 2020

@author: Florian
"""

import numpy as np
from MDPlus.core import *
from scipy.interpolate import interp2d
import mdtraj as md
import sys

def correct_dna(pdb):

    a = open(pdb, 'r')
    f_list = []
    information_residue= []
    information_strand = []
    information_list_pdb = []
    i_residue = []
    i_strand = []
    pdb_list = []     
    pdb_start       = []
    pdb_last        = []

    #Cut dna for obtain dyads and delete some atom of base

    ##only use generic base atoms: limited training set.
    use_atoms = {"O5'", "C5'", "C4'", "O4'", "C1'",\
                 "P", "C5",  "C6",  "N1",  "C2",\
                 "N3", "C4", "C3'", "C2'", "O3'"}

    dont_use_atom_1 = {"C5'", "O5'", "P"}
    dont_use_atom_2 = {"O3'"}
    for line in a:
        as_Liste = line.split(" ")
        as_List = [elem for elem in as_Liste if elem.strip()]
        if as_List[0] == "ATOM":
            #We don't use hydrogen atom for use RX structure or NMR structure or other
            if as_List[2] in use_atoms:
                if as_List[4] in ['A', 'B']: 
                    information_strand.append(as_List[4])
                
                    information_residue.append(int(as_List[len(as_List)-7]))

                    information_list_pdb.append(line)
          
    for i in range(len(information_strand)):
        if information_strand[i] != information_strand[i - 1]:
            i_strand.append(information_strand[i])
    
    #If the strucutre is just one strand, there is zero dna helix
    if len(i_strand) < 2:
        print('zero dna helix')
        sys.exit()
        
    for i in range(len(information_residue)):
        if information_residue[i] != information_residue[i - 1]:
            i_residue.append(information_residue[i])

    #Loop for create a dyads_list for crd format and for pdb format
    for p in range(int(0), int(len(i_residue)/2), int(2)):
        residue_list_pdb = []
        for i in information_list_pdb:
            as_Liste = i.split(" ")
            as_List = [elem for elem in as_Liste if elem.strip()]
            if float(as_List[5]) == i_residue[p] or \
               float(as_List[5]) == i_residue[p + 1]:
                   if as_List[4] == i_strand[0]:
                       residue_list_pdb.append(i)
                       
            if float(as_List[5]) == i_residue[-p - 2] or \
               float(as_List[5]) == i_residue[-p - 1]:
                   if as_List[4] == i_strand[1]: 
                       residue_list_pdb.append(i)
        pdb_list.append(residue_list_pdb)
    np.savetxt('pdb_start.pdb', pdb_list[0], fmt = '%s')
    
    start = md.load('pdb_start.pdb')
    start.save('pdb_start.pdb')
    
    #Delete some other atom in backbone for obtain a dyads with same number of atom. For the first dyads.
    o = open('pdb_start.pdb', 'r')
    for line in o:
        as_Liste = line.split(" ")
        as_List = [elem for elem in as_Liste if elem.strip()]
        if as_List[0] == "ATOM":
            if float(as_List[5]) != i_residue[0] and \
               float(as_List[5]) != i_residue[-2]:
                   if as_List[2] not in dont_use_atom_2:
                       pdb_start.append(line)
                   
            if float(as_List[5]) == i_residue[0] or float(as_List[5]) == i_residue[-2]:
                   if as_List[2] not in dont_use_atom_1:
                       pdb_start.append(line)
                       
    #Sometimes with some dna structure files, some important atom are missing. 
    if len(pdb_start) == 52:
        np.savetxt('dyads_dna[0].pdb', pdb_start, fmt='%s')
        f0 = Fasu('dyads_dna[0].pdb')
        f_list.append(f0)
    
    #Delete some other atom in backbone for obtain a dyads with same number of atom. For all middle dyads.
    for i_dyads_middle in range(1, len(pdb_list) - 1):
        np.savetxt('pdb_middle.pdb', pdb_list[i_dyads_middle], fmt='%s')
        middle = md.load('pdb_middle.pdb')
        middle.save('pdb_middle.pdb')
        o = open('pdb_middle.pdb', 'r')
        pdb_middle = []
        for line in o:
            as_Liste = line.split(" ")
            as_List = [elem for elem in as_Liste if elem.strip()]
            if as_List[0] == "ATOM":
                for i in range(len(i_residue)):
                        if i%2 ==0:
                            if int(as_List[5]) == i_residue[i]:
                                if as_List[2] not in dont_use_atom_1:
                                    pdb_middle.append(line)
                                
                        else:
                            if int(as_List[5]) == i_residue[i]:
                                if as_List[2] not in dont_use_atom_2:
                                    pdb_middle.append(line)
        """                            
        Sometimes with some files this programm take two time each line.
        This loop correct this error.
        If this error doesn't exist, this loop doesn't change the result.
        """
        correct_pdb_middle = []
        for i in range(len(pdb_middle)):
            if pdb_middle[i] != pdb_middle[i - 1]:
                correct_pdb_middle.append(pdb_middle[i])
                
                
        #Sometimes with some dna structure files, some important atom are missing.    
        if len(correct_pdb_middle) == 52:
            np.savetxt('dyads_dna_%f.pdb' % i_dyads_middle, correct_pdb_middle, fmt='%s')
            f_m = Fasu('dyads_dna_%f.pdb' % i_dyads_middle)
            f_list.append(f_m)
    
    np.savetxt('pdb_last.pdb', pdb_list[-1], fmt = '%s')
    
    last = md.load('pdb_last.pdb')
    last.save('pdb_last.pdb')
    
    #Delete some other atom in backbone for obtain a dyads with same number of atom. For the last dyads.
    o = open('pdb_last.pdb', 'r')
    for line in o:
        as_Liste = line.split(" ")
        as_List = [elem for elem in as_Liste if elem.strip()]
        if as_List[0] == "ATOM":
            if int(as_List[5]) != i_residue[int(len(i_residue) / 2) - 2] and \
               int(as_List[5]) != i_residue[int(len(i_residue) / 2)]:
                   if as_List[2] not in dont_use_atom_2:
                       pdb_last.append(line)
                   
            if int(as_List[5]) == i_residue[int(len(i_residue) / 2) - 2] or int(as_List[5]) == i_residue[int(len(i_residue) / 2)]:
                   if as_List[2] not in dont_use_atom_1:
                       pdb_last.append(line)
                       
    #Sometimes with some dna structure files, some important atom are missing.
    if len(pdb_last) == 52:
        np.savetxt('dyads_dna[-1].pdb', pdb_last, fmt='%s')
        f_last = Fasu('dyads_dna[-1].pdb')
        f_list.append(f_last)
    e0_space = np.loadtxt('e0_space_common.txt')
    e1_space = np.loadtxt('e1_space_common.txt')
    correct_E = []
    correct_F = []
    #Calulation of Fasu of reference
    ref_fasu = Fasu('dyads_GC_0.000000.pdb')
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
  
if __name__ == "__main__":
    
    pdb = 'Florian_dna/182D_noInterc.pdb'
    print(correct_dna(pdb))
    

      
    
