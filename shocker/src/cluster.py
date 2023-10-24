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
    clusters the water particles in the system. The resulting clusters are
    collections of bins containing the water particles

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

    def __init__(self, all_atoms, box_dim, bin_size):

        self.all_atoms = all_atoms
        self.box_dim = box_dim
        self.bin_size = bin_size
        self.nr_bins = [int(x/self.bin_size) for x in self.box_dim]

    def bin_converter_l(self, target_lipids):
        """
        Each element of a set of particles with coordinates 'target_lipids' in
        a simulationbox of size 'box_dim' is assigned to a bin of size
        'box_dim/nr_bins'. Units are in Angström.

        Parameters:
        -----------
        target_lipids: array
            coordinates of particles in a simulation box

        Returns:
        --------
        array (x,3) containing the positions of the bins containing the
        particles in 'target_lipids'
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
        Creates a cubic system of bins with dimensions 'nr_bins'. '1' is
        assigned to lipid bins according to 'l_bins'. '0' is assigned to all
        other bins (water)

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
            if i[0] <= self.nr_bins[0]-1 and\
                i[1] <= self.nr_bins[1]-1 and\
                    i[2] <= self.nr_bins[2]-1:

                bsystem[i[0]][i[1]][i[2]] = 1

        return bsystem

    def neighbor_view(self, direction, bin_system, storage, hydration_barrier, cur):
        """
        Of a given bin-system 'system' containing zeros and ones this function
        looks in a particular 'direction' to a 'target' bin relative to the 
        bin 'cur'. If the neighboring bin contains a zero the value is changed 
        to -1. the position of the zero is stored.

        Parameters
        ----------
        direction : tuple (3)
            direction to investigate (for example (1,0,0) is the positive
            x-direction)
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

        target = [cur[0] + x_dir, cur[1] + y_dir, cur[2] + z_dir]

        for i in range(3):
            if target[i] >= self.nr_bins[i]:
                target[i] = target[i] - self.nr_bins[i]
            if target[i] < 0:
                target[i] = target[i] + self.nr_bins[i]

        if bin_system[target[0]][target[1]][target[2]] == 0:
            storage.append(tuple(target))
            bin_system[target[0]][target[1]][target[2]] = -1

        elif bin_system[target[0]][target[1]][target[2]] == 1:
            hydration_barrier.append(tuple(target))

    def cluster_finder(self, bin_system):
        """
        Clusters the zeros in a system 'bin_system 'of zeros and ones. The
        starting point is a zero and looks in 6 directions for another zero,
        which then turns into '-1' and the position is stored. The calculation
        continues until all zeros in this cluster are changed to '-1', after
        which the module looks for a new starting zero elsewehere. Calculation
        is finished as soon as there are no more zeros left in the system.

        Parameters:
        -----------
        bin_system: array (x, y, z, 1)
            cubic system containing zeros and ones

        Returns:
        --------
        array, a list of cluster positions of different sizes
        """
        directions = [(1, 0, 0),
                      (-1, 0, 0),
                      (0, 1, 0),
                      (0, -1, 0),
                      (0, 0, 1),
                      (0, 0, -1)]
        black_list = bin_system
        cluster_list = []
        hb_cluster_list = []

        while zero_finder(black_list) != 0:

            next_start = zero_finder(black_list)

            gray_list = [next_start]
            white_list = []
            temp_list = []
            hb = []
            while len(gray_list) > 0:
                for zero in gray_list:
                    for direction in directions:
                        self.neighbor_view(direction,
                                           black_list,
                                           temp_list,
                                           hb,
                                           zero)

                    white_list.append(zero)
                    black_list[zero[0]][zero[1]][zero[2]] = -1

                gray_list = temp_list
                temp_list = []
            hb_list = white_list + hb
            cluster_list.append(white_list)
            hb_cluster_list.append(hb_list)

        return cluster_list, hb_cluster_list
