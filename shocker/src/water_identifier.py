#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 14:52:13 2022

@author: marco
"""
import numpy as np
import random


class Identifier():
    """
    Selecting water particles to be removed or replaced. these particles are
    chosen from a given water bin cluster

    Parameters:
    -----------
    bin_clusters: array
        a list of bin clusters of the water particles
    box_dim: array
        dimensions of the simulation box
    bin_size: int
        demensions of the bins
    nr_remove: int
        number of removed particles
    water: MDAnalysis atomgroup
        atomgroup containing all solvent particles
    """

    def __init__(self,
                 bin_clusters,
                 box_dim,
                 bin_size,
                 nr_remove,
                 water):

        self.bin_clusters = bin_clusters
        self.box_dim = box_dim
        self.bin_size = bin_size
        self.nr_bins = [int(x/self.bin_size) for x in self.box_dim]
        self.nr_remove = nr_remove
        self.water = water

    def cluster_selecter(self):
        """
        From a list of water clusters the relevant cluster is selected (vesicle
        interior). This is the cluster that does not include bin (0,0,0)
        and the bin in furthest corner from (0,0,0).

        Returns:
        --------
        array (x,3), a list of positions of the desired cluster of bins
        """
        inner_cluster = []
        outer_cluster = []

        for i, _ in enumerate(self.bin_clusters):
            if (0, 0, 0) not in self.bin_clusters[i] and \
                (self.nr_bins[0]-1,
                 self.nr_bins[1]-1,
                 self.nr_bins[2]-1) not in self.bin_clusters[i]:

                inner_cluster = inner_cluster + self.bin_clusters[i]
            else:
                outer_cluster = outer_cluster + self.bin_clusters[i]

        return inner_cluster, outer_cluster

    def bin_converter_w(self, target_atoms):
        """
        Each element of a set of particles with coordinates 'target_atoms' in a
        simulationbox of size 'box_dim' is assigned to a bin of size
        'box_dim/nr_bins'. Units are in Angstr??m.

        Parameters:
        -----------
        target_atoms: array
            coordinates of particles in a simulation box

        Returns:
        --------
        array (x,3) containing the positions of the bins containing the
        particles in 'target_atoms'
        """
        x_pos = np.floor(target_atoms[:, 0]/self.bin_size)
        y_pos = np.floor(target_atoms[:, 1]/self.bin_size)
        z_pos = np.floor(target_atoms[:, 2]/self.bin_size)

        binsfloat = np.stack((x_pos, y_pos, z_pos), axis=-1)
        bins = np.int_(binsfloat)

        for i in bins:
            for pos in range(3):
                if i[pos] > self.nr_bins[pos]:
                    i[pos] = i[pos] - self.nr_bins[pos]
                if i[pos] < 0:
                    i[pos] = self.nr_bins[pos] + i[pos]
        return bins

    def index_finder_m(self, w_all, w_cluster):
        """
        Finds the indices of the clustered water bins 'w_cluster' in the
        complete list of water bins 'w_all'. The search continues until the
        desired number of indices 'nr_remove' is found. In other words: it
        finds the indices of the water particles residing in the clustered bins

        Parameters:
        -----------
        w_all: array (x, 3)
            a list in which each water particle coordinate in the
            system is converted to the position of the bin it resides in.
        w_cluster: array (x, 3)
            a list of positions of clustered water-bins (vesicle interior)

        Returns:
        --------
        array, indices of the water particles we want to remove from the
        vesicle interior according to the total water list
        """
        indices = []

        while len(indices) < self.nr_remove:

            i = random.randint(0, len(w_cluster)-1)
            bin_index = np.where((w_all[:, 0] == w_cluster[i][0]) &
                                 (w_all[:, 1] == w_cluster[i][1]) &
                                 (w_all[:, 2] == w_cluster[i][2]))[0]

            if len(bin_index) > 0:
                i2 = random.randint(0, len(bin_index)-1)
                index = bin_index[i2]
                global_index = self.water[index].index
                if global_index not in indices:
                    indices.append(int(global_index))

        return indices[:self.nr_remove]

    def index_finder_m_aa(self, w_all, w_cluster, oxy_pos):
        """
        Finds the indices of the clustered water bins 'w_cluster' in the
        complete list of water bins 'w_all'. The search continues until the
        desired number of indices 'nr_remove' is found. This specific function
        is used in case of all-atom simulations. A water molecule with its
        oxygen atom inside the bin is seen as 'inside'

        Parameters:
        -----------
        w_all: array (x, 3)
            a list in which each water particle coordinate in the
            system is converted to the position of the bin it resides in.
        w_cluster: array (x, 3)
            a list of positions of clustered water-bins (vesicle interior)
        oxy_pos: integer
            position of oxygen atom relative to the hydrogen atoms in the
            trajectory file (0, 1 or 2)

        Returns:
        --------
        array, indices of the water particles we want to remove from the
        vesicle interior according to the total water list
        """
        indices = []

        while len(indices) < self.nr_remove:

            i = random.randint(0, len(w_cluster)-1)
            bin_index = np.where((w_all[:, 0] == w_cluster[i][0]) &
                                 (w_all[:, 1] == w_cluster[i][1]) &
                                 (w_all[:, 2] == w_cluster[i][2]))[0]
            for index in bin_index:
                condition = (index-oxy_pos)/3
                if condition.is_integer():
                    global_index = self.water[index].index
                    indices.append(int(global_index))

        return indices[:self.nr_remove]
