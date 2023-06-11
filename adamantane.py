#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ase.io import read
from ase.visualize import view
from ase.calculators.emt import EMT
from ase.constraints import FixBondLengths
from ase.geometry.analysis import Analysis
import math
import os
import numpy as np

molecule = read("c10h16.xyz")
molecule.calc = EMT()
original_positions = molecule.get_posisions()

c = FixBondLengths(ch_bonds[0])
molecule.set_constraint(c)

output_file_1 = 'Adamantane: CC_bond length vs Energy.txt'
if not os.path.exists(output_file):
        with open (output_file,'w') as f:
                f.write('{0:<25}{1:<25}\n'.format('distance (A)','energy (eV)')
                
            
                
fac_min =  0.9
fac_max = 1.1
step_size = 0.0025
n_steps = math.ceil((fac_max-fac_min)/step_size)
factors = np.linspace(fac_min,fac_max,n_steps)


for factor in factors:
    
    molecule.set_positions(factor*original_positions)
    cc_distance = molecule.get_distance(0,4)
    energy = molecule.get_potential_energy()
    with open (output_file_1,'a') as f:
        f.write('{0:<25}{1:<25}\n'.format(cc_distance, energy))
        
        

molecule = read("c10h16.xyz")
        
output_file_1 = 'Adamantane Analysis.txt'
if not os.path.exists(output_file):
        with open (output_file,'a') as f:
                f.write('Adamantane: CC_bond length vs Energy\n')
                f.write('\n')

ana_1 = Analysis(molecule)
ch_bonds = ana_1.get_bonds('C', 'H', unique = True)
cc_bonds = ana_1.get_bonds('C', 'C', unique = True)
ch_bond_values = ana_1.get_values(ch_bonds)
cc_bond_values = ana_1.get_values(cc_bonds)
with open (output_file,'a') as f:
    f.write('The number of C-H Bonds are {}\n'.format(len(ch_bonds[0])))
    f.write('The C-H Bonds are {}\n'.format(ch_bonds[0]))
    f.write('The C-H Bond Values are {}\n'.format(ch_bond_values))
    f.write('The number of C-C Bonds are {}\n'.format(len(cc_bonds[0])))
    f.write('The C-C Bonds are {}\n'.format(cc_bonds[0]))    
    f.write('The C-C Bond Values are {}\n'.format(cc_bond_values))
    f.write('\n')

molecule.set_positions(0.9*molecule.get_positions())                
ana_2 = Analysis(molecule)
ch_bonds = ana_2.get_bonds('C', 'H', unique = True)
cc_bonds = ana_2.get_bonds('C', 'C', unique = True)
ch_bond_values = ana_2.get_values(ch_bonds)
cc_bond_values = ana_2.get_values(cc_bonds)
with open (output_file,'a') as f:
    f.write('The number of C-H Bonds are {}\n'.format(len(ch_bonds[0])))
    f.write('The C-H Bonds are {}\n'.format(ch_bonds[0]))
    f.write('The C-H Bond Values are {}\n'.format(ch_bond_values))
    f.write('The number of C-C Bonds are {}\n'.format(len(cc_bonds[0])))
    f.write('The C-C Bonds are {}\n'.format(cc_bonds[0]))    
    f.write('The C-C Bond Values are {}\n'.format(cc_bond_values))
