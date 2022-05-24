#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 09:58:19 2022

@author: marco
"""

import numpy as np

def multiples_remover(allbins):
    """
    Removes multiples appearing in a list of bin positions

    Parameters:
    -----------
    allbins: array (x, 3)
        bin positions, possibly containing multiples of the same bin

    Returns:
    --------
    array (x, 3) of all unique bin positions
    """
    bintuples = [tuple(row) for row in allbins]
    binuniques = np.unique(bintuples, axis=0)

    return binuniques

def zero_finder(bin_system):
    """
    Walking through a cubic bin-system containing zeros and ones, searching for
    a zero.

    Parameters:
    -----------
    bin_system: array (x, y, z, 1)
        cubic system containing zeros and ones

    Returns:
    --------
    array (3, ), the position of the first zero-bin found.
    returns '0' if no zeros are found in the system.
    """
    zero = 0
    for x_pos, _ in enumerate(bin_system):
        for y_pos, _ in enumerate(bin_system[x_pos]):
            for z_pos, _ in enumerate(bin_system[x_pos][y_pos]):
                if bin_system[x_pos][y_pos][z_pos] == 0:
                    zero = (x_pos, y_pos, z_pos)
    return zero

class Cluster():
    """
    clusters the water particles in the system. The resulting clusters are collections
    of bins containing the water particles
    
    Parameters:
    -----------
    lipids: MDAnalysis atomgroup
        atomgroup of the lipids in the system
    box_dim: array
        dimensions of the simulation box
    bin_size: int
        dimensions of the bins
    lipid_list: array
        list of lipid names the vesicle is composed of        
    """
    def __init__(self, lipids, box_dim, bin_size, lipid_list):
        
        self.lipids = lipids
        self.box_dim = box_dim
        self.bin_size = bin_size
        self.nr_bins = [int(x/self.bin_size) for x in self.box_dim]
        self.lipid_list = lipid_list
        
    def l_particle_selector(self, library):
        """
        Of each lipid in 'lipidlist' the last three tail particles are selected
        according to the Martini3 library file.

        Returns:
        --------
        string: the atom selector command for isolating the desired tail atoms
        """
        tail_particles = []

        for lipid in self.lipid_list:

            with open(library, 'r') as lipfile:
                lines = lipfile.readlines()
                for i, _ in enumerate(lines):
                    name = lines[i].strip()[1:-1]

                    if lipid == 'CHOL':
                        tail_particles.append('C1 C2')

                    elif name.replace(' ', '') == lipid:

                        lipid_particles = []
                        counter = 1
                        while lines[i+counter] != '\n':
                            lipid_particles.append(lines[i+counter].split()[1])
                            counter += 1

                        tail_length = int(lipid_particles[-1][1])
                        tail_positions = [tail_length - 2, tail_length - 1, tail_length]

                        for atom in lipid_particles:
                            if atom[1] == str(tail_positions[0]) \
                            or atom[1] == str(tail_positions[1]) \
                            or atom[1] == str(tail_positions[2]):
                                tail_particles.append(atom)

        unique_tail_particles = []
        for i in tail_particles:
            if i not in unique_tail_particles:
                unique_tail_particles.append(i)

        selection_command = 'name'
        for i in unique_tail_particles:
            selection_command = selection_command + ' ' + str(i)
            
        target_lipid_atoms = self.lipids.select_atoms(selection_command)

        return target_lipid_atoms
    
    def l_particle_selector_aa(self):
        """
        In case of all-atom simulations the carbon backbone of the lipid tails
        are used for defining the lipid bins. The names of these particles can be 
        specified in the 'lip' flag
        
        Returns:
        --------
        string: the atom selector command for isolating the desired tail atoms
        """    
        selection_command = 'name'
        for i in self.lipid_list:
            selection_command = selection_command + ' ' + i
        
        print(selection_command)
        target_lipid_atoms = self.lipids.select_atoms(selection_command)
        
        return target_lipid_atoms    

    def bin_converter_l(self, target_lipids):
        """
        Each element of a set of particles with coordinates 'target_lipids' in a
        simulationbox of size 'box_dim' is assigned to a bin of size
        'box_dim/nr_bins'. Units are in Angström.
    
        Parameters:
        -----------
        target_lipids: array
            coordinates of particles in a simulation box
        
        Returns:
        --------
        array (x,3) containing the positions of the bins containing the particles in
        'target_lipids'
        """
        x_pos = np.floor(target_lipids[:, 0]/(self.bin_size))
        y_pos = np.floor(target_lipids[:, 1]/(self.bin_size))
        z_pos = np.floor(target_lipids[:, 2]/(self.bin_size))
    
        binsfloat = np.stack((x_pos, y_pos, z_pos), axis=-1)
        bins = np.int_(binsfloat)
    
        for i in bins:
            for pos in range(3):
                if i[pos] > self.nr_bins[pos]:
                    i[pos] = i[pos] - self.nr_bins[pos]
                if i[pos] < 0:
                    i[pos] = self.nr_bins[pos] + i[pos]
        return bins

    def bin_system_maker(self, l_bins):
        """
        Creates a cubic system of bins with dimensions 'nr_bins'. '1' is assigned
        to lipid bins according to 'l_bins'. '0' is assigned to all other bins
        (water)
    
        Parameters:
        -----------
        l_bins: array (x, 3)
            
        Returns:
        --------
        array (nr_bins[0], nr_bins[1], nr_bins[2], 1), a cubic bin system
        containing zeros (water) and ones (lipids)
        """
        bsystem = np.zeros((self.nr_bins[0], self.nr_bins[1], self.nr_bins[2]))
    
        for i in l_bins:
            if i[0] <= self.nr_bins[0]-1 and i[1] <= self.nr_bins[1]-1 and i[2] <= self.nr_bins[2]-1:
                bsystem[i[0]][i[1]][i[2]] = 1
    
        return bsystem

    def neighbor_view(self, direction, bin_system, storage, cur):
        """
        Of a given bin-system 'system' containing zeros and ones this function looks
        in a particular 'direction' relative to a bin 'cur'. If the neighboring bin
        contains a zero the value is changed to -1. the position of the zero is stored.
    
        Parameters
        ----------
        direction : tuple (3)
            direction to investigate (for example (1,0,0) is the positive x-direction)
        bin_system: array (x, y, z, 1)
            cubic system containing zeros and ones
        storage : array (x,3)
            array in which the position of the found zero-bin is stored
        cur : array (3, )
            the position of the bin from which we observe
    
        Returns
        -------
        Zero is changed to -1 and its position is stored
        """
        x_dir = direction[0]
        y_dir = direction[1]
        z_dir = direction[2]
    
        abs_values = [abs(i) for i in direction]
        axis = abs_values.index(1)
        sign = sum(direction)
    
        if (sign<0 and cur[axis]>0) or (sign>0 and cur[axis]<self.nr_bins[axis]-1):
            if bin_system[cur[0] + x_dir][cur[1] + y_dir][cur[2] + z_dir] == 0:
                storage.append((cur[0] + x_dir, cur[1] + y_dir, cur[2] + z_dir))
                bin_system[cur[0] + x_dir][cur[1] + y_dir][cur[2] + z_dir] = -1

    def cluster_finder(self, bin_system):
        """
        Clusters the zeros in a system 'bin_system 'of zeros and ones. The starting
        point is a zero and looks in 6 directions for another zero, which then
        turns into '-1' and the position is stored. The calculation continues until
        all zeros in this cluster are changed to '-1', after which the module
        looks for a new starting zero elsewehere. Calculation is finished as soon
        as there are no more zeros left in the system.
    
        Parameters:
        -----------
        bin_system: array (x, y, z, 1)
            cubic system containing zeros and ones
        
        Returns:
        --------
        array, a list of cluster positions of different sizes
        """
        directions = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
        black_list = bin_system
        cluster_list = []
    
        while zero_finder(black_list) != 0:
    
            next_start = zero_finder(black_list)
    
            gray_list = [next_start]
            white_list = []
            temp_list = []
    
            while len(gray_list) > 0:
                for zero in gray_list:
                    for direction in directions:
                        self.neighbor_view(direction, black_list, temp_list, zero)
    
                    white_list.append(zero)
                    black_list[zero[0]][zero[1]][zero[2]] = -1
    
                gray_list = temp_list
                temp_list = []
            cluster_list.append(white_list)
        
        return cluster_list

    
