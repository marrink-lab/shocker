#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 09:58:19 2022

@author: marco
"""

import numpy as np

def bin_converter(coord, box_dim, nr_bins):
    """
    Each element of a set of particles with coordinates 'coord' in a
    simulationbox of size 'boxdim' is assigned to a bin of size
    'boxdim/nrbins'. Units are in Angström.

    Parameters:
    -----------
    coord: array (x,3)
        coordinates of particles in a simulation box
    boxdim: array (3,)
        dimensions of the simulation box
    nrbins: array (3,)
        number of bins in the x-, y- and z-direction

    Returns:
    --------
    array (x,3) containing the positions of the bins housing the particles in
    'coord'
    """
    x_pos = np.ceil(coord[:, 0]/(box_dim[0]/nr_bins[0]))
    y_pos = np.ceil(coord[:, 1]/(box_dim[1]/nr_bins[1]))
    z_pos = np.ceil(coord[:, 2]/(box_dim[2]/nr_bins[2]))

    binsfloat = np.stack((x_pos, y_pos, z_pos), axis=-1)
    bins = np.int_(binsfloat)

    for i in bins:
        for pos in range(3):
            if i[pos] > nr_bins[pos]:
                i[pos] = i[pos] - nr_bins[pos]
            if i[pos] < 0:
                i[pos] = nr_bins[pos] + i[pos]
    return bins

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

def bin_system_maker(l_bins, nr_bins):
    """
    Creates a cubic system of bins with dimensions 'nr_bins'. '1' is assigned
    to lipid bins according to 'l_bins'. '0' is assigned to all other bins
    (water)

    Parameters:
    -----------
    l_bins: array (x, 3)
        the positions of lipid bins
    nr_bins: array (3, )
        number of bins in the x-, y- and z-direction

    Returns:
    --------
    array (nr_bins[0], nr_bins[1], nr_bins[2], 1), a cubic bin system
    containing zeros (water) and ones (lipids)
    """
    bsystem = np.zeros((nr_bins[0], nr_bins[1], nr_bins[2]))

    for i in l_bins:
        bsystem[i[0]-1][i[1]-1][i[2]-1] = 1

    return bsystem

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

def index_finder(w_all, w_cluster, nr_remove):
    """
    Finds the indices of the clustered water bins 'w_cluster' in the complete
    list of water bins 'w_all'. The search continues until the desired number
    of indices 'nr_remove' is found. In other words: it finds the indices of the
    water particles residing in the clustered bins

    Parameters:
    -----------
    w_all: array (x, 3)
        a list in which each water particle coordinate in the
        system is converted to the position of the bin it resides in.
    w_cluster: array (x, 3)
        a list of positions of clustered water-bins (vesicle interior)
    nr_remove: integer
        the number of water particles we want to remove from
        the vesicle interior

    Returns:
    --------
    array (x, 3), indices of the water particles we want to remove from the
    vesicle interior according to the total water list
    """
    indices = []
    i = 0
    while len(indices) < nr_remove:
        bin_index = np.where((w_all[:, 0] == w_cluster[i][0])\
                             & (w_all[:, 1] == w_cluster[i][1])\
                                 & (w_all[:, 2] == w_cluster[i][2]))[0]
        for index in bin_index:
            indices.append(int(index))
        i = i + 1

    return indices[:nr_remove]

def neighbor_view(nr_bins, direction, bin_system, storage, cur):
    """
    Of a given bin-system 'system' containing zeros and ones this function looks
    in a particular 'direction' relative to a bin 'cur'. If the neighboring bin
    contains a zero the value is changed to -1. the position of the zero is stored.

    Parameters
    ----------
    nr_bins : array (3, )
        number of bins in the x-, y- and z-direction
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
    Zero is changed to one and its position is stored
    """
    x_dir = direction[0]
    y_dir = direction[1]
    z_dir = direction[2]

    abs_values = [abs(i) for i in direction]
    axis = abs_values.index(1)
    sign = sum(direction)

    if (sign<0 and cur[axis]>0) or (sign>0 and cur[axis]<nr_bins[axis]-1):
        if bin_system[cur[0] + x_dir][cur[1] + y_dir][cur[2] + z_dir] == 0:
            storage.append((cur[0] + x_dir, cur[1] + y_dir, cur[2] + z_dir))
            bin_system[cur[0] + x_dir][cur[1] + y_dir][cur[2] + z_dir] = -1

def cluster_finder(bin_system, nr_bins):
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
    nr_bins: array (3, )
        number of bins in the x-, y- and z-direction

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
                    neighbor_view(nr_bins, direction, black_list, temp_list, zero)

                white_list.append(zero)
                black_list[zero[0]][zero[1]][zero[2]] = -1

            gray_list = temp_list
            temp_list = []
        cluster_list.append(white_list)


    return cluster_list

def l_particle_selector(lipidlist):
    """
    Of each lipid in 'lipidlist' the last three tail particles are selected
    according to the Martini library file.

    Parameters:
    -----------
    lipidlist: array
        list of lipid names the vesicle membrane is composed of.

    Returns:
    --------
    string: the atom selector command for isolating the desired tail atoms
    """
    tail_particles = []

    for lipid in lipidlist:

        with open('Martini3.LIB', 'r') as lipfile:
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

    return selection_command
