#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:55:32 2022

@author: marco
"""

class Remover():
    
    """
    Removing water particles from the system
    
    Parameters:
    -----------
    top_old: string
        name of the old topology file
    all_atoms: MDAnalysis atomgroup
        atomgroup of all the particles in the system
    nr_removed: int
        number of removed particles
    """
    
    def __init__(self, top_old, all_atoms, nr_removed, gro_file):
        
        self.top_old = top_old
        self.all_atoms = all_atoms
        self.nr_removed = nr_removed
        self.gro_file = gro_file

    def water_remover_gro(self, indices):
        """
        A string is created of particle indices that have to be removed from the
        system, according to which the new gro file is constructed.
    
        Parameters:
        -----------
        indices: array 
            indices of the particles to be removed
        gro_file: string
            general name of the gro file fed to the system
    
        Returns:
        --------
        writes a new gro file
        """
        nstring = 'not index '
        for i in indices:
            nstring = nstring + str(i) + ' '
        
        keep_atoms = self.all_atoms.select_atoms(nstring)
        keep_atoms.write(self.gro_file)
               
    def water_remover_top(self, new):
        """
        Creates a new topology file in which the number of water particles is modified.
    
        Parameters:
        -----------
        new: string
            path to the new topology file
           
        Returns:
        --------
        writes new topology file
        """
        with open(self.top_old) as cur_top:
            lines_top = cur_top.readlines()
            nr_top_lines = len(lines_top)
    
    
            for i in range(nr_top_lines):
                sp_string_top = lines_top[i].split(' ')
                if sp_string_top[0] != 'W':
                    with open(new, 'a') as new_top:
                        new_top.write(lines_top[i])
                else:
    
    
                    cur_nr_water = sp_string_top[1].split('\n')
                    new_nr_water = int(cur_nr_water[0]) - self.nr_removed
                    new_string_element_water = str(new_nr_water) + '\n'
    
                    sp_string_top[1] = new_string_element_water
                    new_string_water = ' '.join(sp_string_top)
    
                    with open(new, 'a') as new_top:
                        new_top.write(new_string_water)