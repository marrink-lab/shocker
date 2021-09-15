#!/usr/bin/env python
# coding: utf-8

# In[31]:


import MDAnalysis as mda
import numpy as np
import subprocess

  


# In[34]:


def water_selecter():
        
    # Calculating the center of mass (z-coordinate) of both bilayers
    top_atoms = lipid_atoms[int(nr_lipid_atoms/2):int(nr_lipid_atoms)]
    bottom_atoms = lipid_atoms[0:int(nr_lipid_atoms/2)]

    upper_limit_z = top_atoms.centroid()[2]
    lower_limit_z = bottom_atoms.centroid()[2]

    # Finding the indices of the water particles situated between the two bilayers
    water_inside = []

    for i in range(nr_water_atoms):
        water_z_pos = water_atoms.positions[i][2]
        if lower_limit_z < water_z_pos < upper_limit_z:
            water_inside.append(i)
       
    # Selecting the atoms to be removed
    removed_atoms = water_inside[0:nr_removed_atoms]             # Local indices of removed atoms
    print(len(water_inside))
    
    return removed_atoms
    


# In[28]:


def water_remover_gro(removed):
    
    removed_atoms_li = removed
    
    # The following piece of code copies the original .gro file in which the desired percentage of water particles 
    # is removed.

    removed_atoms_gi = []                              # Global indices of removed atoms
    for i in removed_atoms_li:
        removed_atoms_gi.append(i + nr_lipid_atoms)
    removed_atoms_gi.append(-1)                        # -1 has to be added to make sure the list is never empty

    cur_gro = open(original_gro)
    gro_lines = cur_gro.readlines()
    new_gro = open(modified_gro,'a')

    # New .gro file is created in 4 steps:

    # Step 1: The first line (title)
    new_gro.write(gro_lines[0])

    # Step 2: The total number of particles in the system. The number of removed particles has to be subtracted
    split_number = gro_lines[1].split('\n')
    new_number = int(split_number[0]) - nr_removed_atoms
    new_string = str(new_number) + '\n'
    new_gro.write(new_string)

    # Step 3: All lines (particles) of the original file are copied except for the particles that have to be removed.
    # Every copied particle gets a new index number to yield a continuous list of indexes in the file.
    index = 1
    nr_lines = len(gro_lines)
    for i in range(2, nr_lines - 1):
        if i != removed_atoms_gi[0] + 2:                  # +2 is needed since we start at the third line
            sp_string = gro_lines[i].split(' ')           # The line is split to separate the elements the line is
            cur_pos = len(sp_string) - 1                  # composed of. We start at the end of the line and search for
            counter = 7                                   # the position 7th element (the index number), after which
            while counter != 0:                           # this number is changed into the current 'index' value.
                if sp_string[cur_pos] != '':
                    counter = counter - 1
                cur_pos = cur_pos - 1
            index_pos = cur_pos + 1

            index_length = len(str(index))                     # We need to know the number of character of the current index
            if index_length == 5:                              # value, because at a length of 5 characters we need a different
                atom_name = sp_string[index_pos][0:-5]         # approach. In this case the atom name and the index number
                sp_string[index_pos] = atom_name + str(index)  # are 1 string, which has to be seperated first, and rejoined.

            else:
                sp_string[index_pos] = str(index)

            new_line = ' '.join(sp_string)
            new_gro.write(new_line)
            index = index + 1

        else:
            removed_atoms_gi.pop(0)

    # Step 4: The last line with the current dimensions of the box is added.       
    new_gro.write(gro_lines[-1])

    cur_gro.close()
    new_gro.close()


# In[7]:


def water_remover_top():
    
    # Next, the original topology file is modified to make sure it shows the correct number of W-particles
    cur_top = open(original_top)
    lines_top = cur_top.readlines()
    new_top = open(modified_top, 'a')
    nr_top_lines = len(lines_top)

    # Step 1, all lines are copied to the new .top file except for the last one, containing the number of W-particles.
    for i in range(nr_top_lines - 1):
        new_top.write(lines_top[i])

    # Step 2, the number of removed W-particles has to be subtracted from the current number.
    sp_string_top = lines_top[-1].split(' ')
    cur_nr_water = sp_string_top[1].split('\n')
    new_nr_water = int(cur_nr_water[0]) - nr_removed_atoms
    new_string_element_water = str(new_nr_water) + '\n'

    sp_string_top[1] = new_string_element_water
    new_string_water = ' '.join(sp_string_top)

    new_top.write(new_string_water)

    new_top.close()


# In[ ]:


shock = 0
while shock != 16:
    subprocess.call('gmx grompp -f martini_md.mdp -c POPC_minimized_eq.gro -p topol2.top -o POPC_minimized_eq2.tpr', shell=True)
    subprocess.call('gmx mdrun -s POPC_minimized_eq2.tpr -v -x POPC_minimized_eq2.xtc -c POPC_minimized_eq2.gro', shell=True)
    subprocess.call('rm POPC_minimized_eq.gro', shell=True)
    
    original_gro = 'POPC_minimized_eq2.gro'
    modified_gro = 'POPC_minimized_eq.gro'
    original_top = 'topol2.top'
    modified_top = 'topol.top'
    u = mda.Universe(original_gro)
    lipid_atoms = u.select_atoms('resname POPC')
    water_atoms = u.select_atoms('resname W')
    nr_lipid_atoms = len(lipid_atoms)
    nr_water_atoms = len (water_atoms)
    nr_removed_atoms = 400
    
    removed = water_selecter()
    water_remover_gro(removed)
    water_remover_top()
    subprocess.call('rm topol2.top', shell=True)
    subprocess.call('mv topol.top topol2.top', shell=True)
    shock = shock + 1 

