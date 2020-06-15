# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:35:20 2020

@author: Florian
"""

import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj as md


pdb = 'Florian_dna/ggcgac_plainFix_dnaOnly.pdb'
top = 'Florian_dna/ggcgac_plainFix_dnaOnly.top'
nc = 'Florian_dna/ggcgac_plainFix_dnaOnly.nc'

def correction_DNA_energy(pdb, top, nc):

    from correct_DNA import correct_dna    
    correct = []
    #Load the trajectory
    traj = md.load(nc, top = pdb)
    #Selection of each frame of trajectory
    for i in range(len(traj)):
        for ii in range(len(traj[i])): 
            #Creation of files for calculation for each frame
            traj[i][ii].save('frame_DNA.ncrst')
            traj[i][ii].save('frame_DNA.pdb')
            prmtop = AmberPrmtopFile(top)
            inpcrd = AmberInpcrdFile('frame_DNA.ncrst')
            ##we shouldn't need an integrator to do a single-point energy evaluation
            integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

            system = prmtop.createSystem(nonbondedMethod=NoCutoff,
                                     constraints=None)
            simulation = Simulation(prmtop.topology, system, integrator)
            simulation.context.setPositions(inpcrd.positions)
            state      = simulation.context.getState(getEnergy=True, getForces=True)
            ene = state.getPotentialEnergy()
            ene = ene/(kilocalorie/mole)
            #Calculation of energy correction
            correction = correct_dna('frame_DNA.pdb')
            #Creation of correct energy list
            correct.append(ene + correction[0])
    return correct

print(correction_DNA_energy(pdb, top, nc))