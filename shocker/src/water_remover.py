#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 14:52:13 2022

@author: marco
"""
import subprocess
import numpy as np
import Shocker.shocker.src.cluster as cl

def cluster_selecter(clusterlist, nr_bins):
    """
    From a list of water clusters the relevant cluster is selected (vesicle
    interior). This is the cluster that does not include bin (0,0,0) and the
    bin in furthest corner from (0,0,0).

    Parameters:
    -----------
    clusterlist: array
        a list of cluster positions of different sizes
    nr_bins: array (3,)
        number of bins in the x-, y- and z-direction

    Returns:
    --------
    array (x,3), a list of positions of the desired cluster of bins
    """
    target_cluster = []
    for i, _ in enumerate(clusterlist):
        if (0, 0, 0) not in clusterlist[i] and \
            (nr_bins[0]-1, nr_bins[1]-1, nr_bins[2]-1) not in clusterlist[i]:
            target_cluster = target_cluster + clusterlist[i]

    return target_cluster

def water_selecter(lipids, water, nr_bins, box_dim, nr_remove):
    """
    A cubic bin-system in which a 1 is assigned to a lipid bin and a 0 is
    assigned to a water bin. This system is clustered, after which a part of
    the clustered water bins of the vesicle interior are translated back to
    real particle coordinates, which serve as particles that we want to remove
    from the vesicle.

    Parameters:
    -----------
    lipids: array (x,3)
        coordinates of the lipid atoms that have to be targeted
        when assigning lipid bins
    water: array (x,3)
        the positions of the water particles
    nr_bins: array (3,)
        number of bins in the x-, y- and z-direction
    box_dim: array (3,)
        dimensions of the simulation box
    nr_remove: integer
        number of water particles that have to be removed from
        the vesicle interior

    Returns:
    --------
    array, list of indices of water particles that will be removed from the
        vesicle. These indices are local, meaning that these are relative to the
        water particle list only.
    """
    lipid_bins = cl.bin_converter(lipids, box_dim, nr_bins)

    lipid_bins_single = cl.multiples_remover(lipid_bins)

    bin_system = cl.bin_system_maker(lipid_bins_single, nr_bins)

    clusters = cl.cluster_finder(bin_system, nr_bins)

    inner_water_cluster = cluster_selecter(clusters, nr_bins)

    np.random.shuffle(inner_water_cluster)

    water_bins = cl.bin_converter(water, box_dim, nr_bins)

    remove = cl.index_finder(water_bins, inner_water_cluster, nr_remove)

    return remove

def local_to_global(no_water, local_index):
    """
    Local indices of the water particles that have to be removed are converted
    to global indices of the all atoms list.

    Parameters:
    -----------
    no_water: integer
        the number of particles other than water
    local_index: array
        the local indices of water particles to be removed from
        the vesicle

    Returns:
    --------
    array, the global indices of water particles to be removed from the vesicle
    """
    global_index = []
    for i in local_index:
        global_index.append(i + no_water)

    return global_index

def not_string(positions):
    """
    Produces a string used in mdanalysis to construct the new gro file without
    the removed water particles.

    Parameters:
    ----------
    positions: array
        A list of indices of the to be removed water particles.

    Returns:
    -------
    nstring : string
        string to use in mdanalysis to exclude the removed water particles
    """
    nstring = 'not index '
    for i in positions:
        nstring = nstring + str(i) + ' '

    return nstring

def water_remover_gro(positions, universe, gro_file):
    """
    A string is created of particle indices that have to be removed from the
    system, according to which the new gro file is constructed.

    Parameters:
    -----------
    positions: array
        indices of atoms to be removed
    universe:
        contains all the information describing the molecular dynamics system.
    gro_file: string
        general name of the gro file fed to the system

    Returns:
    --------
    writes a new gro file
    """
    nstring = 'not index '
    for i in positions:
        nstring = nstring + str(i) + ' '

    keep_atoms = universe.select_atoms(nstring)
    keep_atoms.write(gro_file)

def water_remover_top(old, new, nr_removed):
    """
    Creates a new topology file in which the number of water particles is modified.

    Parameters:
    -----------
    old: string
        path to the old topology file
    new: string
        path to the new topology file
    nr_removed: integer
        number of removed water particles

    Returns:
    --------
    writes new topology file
    """
    with open(old) as cur_top:
        lines_top = cur_top.readlines()
        nr_top_lines = len(lines_top)


        for i in range(nr_top_lines):
            sp_string_top = lines_top[i].split(' ')
            if sp_string_top[0] != 'W':
                with open(new, 'a') as new_top:
                    new_top.write(lines_top[i])
            else:


                cur_nr_water = sp_string_top[1].split('\n')
                new_nr_water = int(cur_nr_water[0]) - nr_removed
                new_string_element_water = str(new_nr_water) + '\n'

                sp_string_top[1] = new_string_element_water
                new_string_water = ' '.join(sp_string_top)

                with open(new, 'a') as new_top:
                    new_top.write(new_string_water)

def name_generator(old, shock_nr):
    """
    Generates a unique name for .gro,.top and .tpr files according to the cycle
    number (shock_nr).

    Parameters:
    -----------
    old: string
        old file name
    shock_nr: integer
        shock number

    Returns:
    --------
    file is saved under the new name in folder 'shockfiles'
    """
    no_ext = old.split('.')[0]
    ext = old.split('.')[1]
    new_name = no_ext + '_' + 's' + str(shock_nr) + '.' + ext

    name_change_command = 'mv ' + old + ' ' + new_name
    subprocess.call(name_change_command, shell=True)
    dir_change = 'mv ' + new_name + ' shockfiles'
    subprocess.call(dir_change, shell=True)
