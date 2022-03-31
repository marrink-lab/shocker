#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 10:04:30 2022

@author: marco
"""

import MDAnalysis as mda
import random
import numpy as np 
from math import dist

def index_finder(w_all, w_bin):
    """
    Finds the indices of the particles in a single bin 'w_bin' in the complete
    list of water bins 'w_all'.

    Parameters:
    -----------
    w_all: array (x, 3)
        a list in which each particle coordinate in the
        system is converted to the position of the bin it resides in.
    w_bin: array (x, 3)
        a list of positions of clustered water-bins

    Returns:
    --------
    array (x, 3), indices of the water particles we want to remove from the
    vesicle interior according to the total water list
    """
    indices = []

    bin_index = np.where((w_all[:, 0] == w_bin[0])\
                         & (w_all[:, 1] == w_bin[1])\
                             & (w_all[:, 2] == w_bin[2]))[0]
    for index in bin_index:
        indices.append(int(index))


    return indices

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
    def __init__(self, outer_cluster, all_atoms, bin_size, nr_removed, box_dim):
        self.outer_cluster = outer_cluster
        self.all_atoms = all_atoms
        self.bin_size = bin_size
        self.box_dim = box_dim
        self.nr_bins = [int(x/self.bin_size) for x in self.box_dim]
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

    def position_generator(self, chosen_bins, all_bins):
        """
        in each bin in "chosen_bins" the 'best' position for a particle is found
        by randomly picking a position and calculating the shortest distance between
        this position and the other particles in the bin. this process is repeated
        10.000 times.

        Parameters:
        -----------
        chosen_bins: array
            list of bin positions
        all_bins: array
            all particles in the system converted to bins

        Returns:
        --------
        array containing positions for water particles to be replaced
        """
        placement_pos = []
        for b in chosen_bins:

            indices = index_finder(all_bins, b)
            istring = 'index '
            for i in indices:
                istring = istring + str(i) + ' '

            w_particles = self.all_atoms.select_atoms(istring).positions
            cur_min = 0
            best_pos = 0
            c = 0
            while c < 10000:
                randx = random.uniform((b[0]+0.2)*self.bin_size, (b[0]+0.8)*self.bin_size)
                randy = random.uniform((b[1]+0.2)*self.bin_size, (b[1]+0.8)*self.bin_size)
                randz = random.uniform((b[2]+0.2)*self.bin_size, (b[2]+0.8)*self.bin_size)
                temp_pos = [randx, randy, randz]

                distance = []
                for particle in w_particles:
                    distance.append(dist(particle, temp_pos))
                if min(distance) > cur_min:
                    cur_min = min(distance)
                    best_pos = temp_pos
                c = c + 1
            print(cur_min)

            placement_pos.append(best_pos)

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