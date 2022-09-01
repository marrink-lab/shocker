#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 10:04:30 2022

@author: marco
"""

import random
import numpy as np
import numpy.linalg


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

    bin_index = np.where((w_all[:, 0] == w_bin[0]) &
                         (w_all[:, 1] == w_bin[1]) &
                         (w_all[:, 2] == w_bin[2]))[0]
    for index in bin_index:
        indices.append(int(index))

    return indices


class Mover():
    """
    Selecting bins in which a random position is chosen to replace water
    particles

    Parameters:
    -----------
    receiving_cluster: array
        locations of the bins of the cluster water particles are moved into
    all_atoms: MDAnalysis atomgroup
        atomgroup of all the particles in the system
    bin_size: int
        dimensions of the bins
    nr_removed: int
        number of removed particles
    box_dim: array
        dimensions of the simulation box
    all_bins: array
        each particle position translated to bin location
    water_group: MDAnalysis atomgroup
        all particles belonging to the solvent
    """

    def __init__(self,
                 receiving_cluster,
                 all_atoms,
                 bin_size,
                 nr_removed,
                 box_dim,
                 all_bins,
                 water_group):

        self.receiving_cluster = receiving_cluster
        self.all_atoms = all_atoms
        self.bin_size = bin_size
        self.box_dim = box_dim
        self.nr_bins = [int(x/self.bin_size) for x in self.box_dim]
        self.nr_removed = nr_removed
        self.all_bins = all_bins
        self.water_group = water_group

    def repositioning_bin_identifier(self):
        """
        Randomly select bins to replace water particles. each water particle
        ends up in a unique bin that only contains water particles. If no bins
        are found containing solely water within reasonable time, the best
        alternative option is used.

        Returns:
        --------
        array containing the bins water particles are moved to
        """
        cluster_length = len(self.receiving_cluster)

        chosen_bins = []
        rand_nr_storage = []
        best_not_w_count = 0
        best_bin = 0
        tries = 0
        while len(chosen_bins) < self.nr_removed:
            random_nr = random.randint(0, cluster_length-1)

            if random_nr not in rand_nr_storage:

                rand_nr_storage.append(random_nr)
                chosen_bin = self.receiving_cluster[random_nr]

                not_w_count = 0

                indices = index_finder(self.all_bins, chosen_bin)

                for i in indices:
                    istring = 'index ' + str(i)
                    if len(self.water_group.select_atoms(istring)) == 0:
                        not_w_count = not_w_count + 1

                if not_w_count == 0:
                    chosen_bins.append(chosen_bin)
                    best_not_w_count = 0
                else:
                    if not_w_count < best_not_w_count:
                        best_not_w_count = not_w_count
                        best_bin = self.receiving_cluster[random_nr]
                        tries = tries + 1
                if tries == 10:
                    chosen_bins.append(best_bin)
                    tries = 0

        return chosen_bins

    def position_generator(self, chosen_bins):
        """
        in each bin in "chosen_bins" the 'best' position for a particle is
        found by randomly picking a position and calculating the shortest
        distance between this position and the other particles in the bin.
        This process is repeated 10.000 times.

        Parameters:
        -----------
        chosen_bins: array
            list of bin positions

        Returns:
        --------
        array containing positions for water particles to be replaced
        """
        sum_min_dist = 0
        placement_pos = []
        for b in chosen_bins:

            indices = index_finder(self.all_bins, b)
            istring = 'index '
            for i in indices:
                istring = istring + str(i) + ' '

            w_particles = self.all_atoms.select_atoms(istring).positions

            cur_min = 0
            best_pos = 0
            c = 0
            while c < 10000:
                randx = random.uniform((b[0]+0.2)*self.bin_size,
                                       (b[0]+0.8)*self.bin_size)
                randy = random.uniform((b[1]+0.2)*self.bin_size,
                                       (b[1]+0.8)*self.bin_size)
                randz = random.uniform((b[2]+0.2)*self.bin_size,
                                       (b[2]+0.8)*self.bin_size)
                temp_pos = [randx, randy, randz]

                distance = []
                for particle in w_particles:
                    distance.append(np.linalg.norm(particle - temp_pos))
                if min(distance) > cur_min:
                    cur_min = min(distance)
                    best_pos = temp_pos
                c = c + 1

            sum_min_dist = sum_min_dist + cur_min

            placement_pos.append(best_pos)

        mean_min_dist = sum_min_dist/len(placement_pos)

        return placement_pos, mean_min_dist

    def position_generator_aa(self, chosen_bins, oxy_pos):
        """
        in each bin in "chosen_bins" the 'best' position for a particle is
        found by randomly picking a position and calculating the shortest
        distance between this position and the other particles in the bin.
        This process is repeated 10.000 times. This function is used in
        all-atom simulations where the center of geometry of water molecules
        is used to relocate them.

        Parameters:
        -----------
        chosen_bins: array
            list of bin positions
        oxy_pos: integer
            position of oxygen atom relative to the hydrogen atoms in the
            trajectory file (0, 1 or 2)

        Returns:
        --------
        array containing positions for water particles to be relocated
        """
        sum_min_dist = 0
        placement_pos = []
        for b in chosen_bins:

            indices = index_finder(self.all_bins, b)
            istring = 'index '
            for i in indices:
                istring = istring + str(i) + ' '

            ox_in_bin = self.all_atoms.select_atoms(istring).positions

            cur_min = 0
            best_pos = 0
            c = 0
            while c < 1000:
                randx = random.uniform((b[0]+0.2)*self.bin_size,
                                       (b[0]+0.8)*self.bin_size)
                randy = random.uniform((b[1]+0.2)*self.bin_size,
                                       (b[1]+0.8)*self.bin_size)
                randz = random.uniform((b[2]+0.2)*self.bin_size,
                                       (b[2]+0.8)*self.bin_size)
                temp_pos = [randx, randy, randz]

                distance = []
                for particle in ox_in_bin:
                    distance.append(dist(particle, temp_pos))
                if min(distance) > cur_min:
                    cur_min = min(distance)
                    best_pos = temp_pos
                c = c + 1

            sum_min_dist = sum_min_dist + cur_min
            placement_pos.append(best_pos)

        mean_min_dist = sum_min_dist/len(placement_pos)

        return placement_pos, mean_min_dist

    def water_repositioning_gro(self, indices, new_pos, gro_file):
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

    def water_repositioning_gro_aa(self, indices, new_pos, gro_file, oxy_pos):
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
        oxy_pos: integer
            position of oxygen atom relative to the hydrogen atoms in the
            trajectory file (0, 1 or 2)

        Returns:
        --------
        new gro file
        """
        pos = self.all_atoms.positions

        for i, e in enumerate(indices):
            if oxy_pos == 0:
                istring = 'index {}'.format(e)
                ox_water = self.all_atoms.select_atoms(istring).\
                    center_of_geometry()
                diff = new_pos[i] - ox_water
                pos[e] = pos[e] + diff
                pos[e+1] = pos[e+1] + diff
                pos[e+2] = pos[e+2] + diff

            if oxy_pos == 1:
                istring = 'index {} {} {}'.format(e, e+1, e-1)
                cog_water = self.all_atoms.select_atoms(istring).\
                    center_of_geometry()
                diff = new_pos[i] - cog_water
                pos[e] = pos[e] + diff
                pos[e+1] = pos[e+1] + diff
                pos[e-1] = pos[e-1] + diff

            if oxy_pos == 2:
                istring = 'index {} {} {}'.format(e, e-1, e-2)
                cog_water = self.all_atoms.select_atoms(istring).\
                    center_of_geometry()
                diff = new_pos[i] - cog_water
                pos[e] = pos[e] + diff
                pos[e-1] = pos[e-1] + diff
                pos[e-2] = pos[e-2] + diff

        self.all_atoms.positions = pos
        self.all_atoms.write(gro_file)
