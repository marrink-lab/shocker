#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 10:04:30 2022

@author: marco
"""

import MDAnalysis as mda
import random 

class Mover():
    """
    Selecting bins in which a random position is chosen to replace water particles
    
    Parameters:
    -----------
    outer_cluster: array
        locations of the bins of the outer water cluster
    all_atoms: MDAnalysis atomgroup
        atomgroup of all the particles in the system
    bin_size: int
        dimensions of the bins
    nr_removed: int
        number of removed particles
    """
    def __init__(self, outer_cluster, all_atoms, bin_size, nr_removed):
        self.outer_cluster = outer_cluster
        self.all_atoms = all_atoms
        self.bin_size = bin_size
        self.nr_removed = nr_removed

    def replacement_bin_identifier(self):
        """
        Randomly select bins to replace water particles
        """
        cluster_length = len(self.outer_cluster)
        
        chosen_bins_nr = []
        for i in range(self.nr_removed):
            random_nr = random.randint(0, cluster_length-1)
            chosen_bins_nr.append(random_nr)
        
        chosen_bins = []
        for i in chosen_bins_nr:
            chosen_bins.append(self.outer_cluster[i])
        
        return chosen_bins
    
    def position_generator(self, bins):
        """
        in each bin in "bins" a position is randomly picked to place one water particle
        
        Parameters:
        -----------
        bins: array
            list of bin positions

        Returns:
        --------
        array containing positions for water particles to be replaced
        """
        placement_pos = []
        for i in bins:
            
            randx = random.uniform((i[0]-1)*self.bin_size, i[0]*self.bin_size)
            randy = random.uniform((i[1]-1)*self.bin_size, i[1]*self.bin_size)
            randz = random.uniform((i[2]-1)*self.bin_size, i[2]*self.bin_size)
            placement_pos.append([randx,randy,randz])
            
        return placement_pos
    
    def water_replacement_gro(self, indices, new_pos, gro_file):
        """
        new positions are assigned to water particles that have to be moved, 
        subsequently a new gro file is created
        
        Parameters:
        -----------
        indices: array
            the indices of the particles to be replaced
        new_pos: array
            the new positions of the particles to be replaced
        gro_file: string
            name of the new gro file
            
        Returns:
        --------
        new gro file
        """
        pos = self.all_atoms.positions
        for i, e in enumerate(indices):
            pos[e] = new_pos[i]
        
        self.all_atoms.positions = pos
        self.all_atoms.write(gro_file)
    

        
        
