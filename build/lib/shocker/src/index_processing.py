# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 11:05:18 2022

@author: s103567
"""

import numpy as np

def ndx_to_ag(all_atoms, infile):
    """
    Create MDAnalysis AtomGroups from a Universe and an index file.
    Examples
    --------
    .. code-block::
    u = mda.Universe(TPR, XTC)
    with open('index.ndx') as infile:
        for group_name, atom_group in ndx_to_ag(u, infile):
            print(group_name, atom_group)
    Parameters
    ----------
    infile: file
        The open index file to read.
    Yields
    ------
    group_name: str
        The name of the group as written in the index file.
    atoms: mda.AtomGroup
        An atom group containing the atoms for that group.
    """
    atoms = all_atoms.atoms
    group_name = None
    indices = []
    for line in infile:
        comment_idx = line.find(';')
        if comment_idx >= 0:
            line = line[comment_idx:]
        line = line.strip()
        if line.startswith('['):
            if group_name is not None:
                indices = np.array(indices, dtype=int) - 1
                yield (group_name, atoms[indices])
            group_name = line[1:-1].strip()
            indices = []
        else:
            indices += line.split()
    if group_name is not None:
        indices = np.array(indices, dtype=int) - 1
        yield (group_name, atoms[indices])
        
def ag_extractor(all_atoms, index_file):
    """
    from the atomgroups available in the index file the coordinates of 
    the membrane particles or atoms are extracted
    
    Returns
    -------
    mem_pos : numpy array
        positions of the membrane particles.

    """
    with open(index_file) as infile:
        agroup = []
        name = []
        for k, l in ndx_to_ag(all_atoms, infile):
            agroup.append(l)
            name.append(k)
             
    return agroup, name

def ag_finder(ag_name, ag_list):
    """
    Returns the atomgroup that matches the given name
    """
    index = ag_list[1].index(ag_name)
    ag = ag_list[0][index]
                   
    return ag

def name_concatenator(name_list):
    """
    concatenates a list of atomgroup names, used to define the groups in the
    mdp files.
    
    Parameters:
    -----------
    name_list: array
    
    Returns:
    --------
    string with group names
    """
    name_str = name_list[0]
    for i in range(1, len(name_list)):
        name_str = name_str + ' ' + name_list[i]
        
    return name_str  