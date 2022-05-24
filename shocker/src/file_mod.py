#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 18:02:38 2022

@author: marco
"""
import subprocess

def name_generator(old, shock_nr, test='no'):
    """
    Generates a unique name for .gro, .top and .tpr files according to the cycle
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
    
    if test == 'no':
        name_change_command = 'mv ' + old + ' ' + new_name
        subprocess.call(name_change_command, shell=True)
        dir_change = 'mv ' + new_name + ' shockfiles'
        subprocess.call(dir_change, shell=True)
    else:
        return new_name
    
def mdp_value_changer(old, new, parameter, value, test='no'):
    """
    Modifies a value in the mdp file by writing a new file and deleting the old
    one.
    
    Parameters:
    -----------
    old: string
        path to the old mdp file
    new: string
        path to the new mdp file
    parameter: string
        parameter name to be modified
    value: integer
        new parameter value
    test: string
        name change and file replacement not exectuted if a test is run
        
    Returns:
    --------
    writes new mdp file with modified parameters
    """
    with open(old) as cur_mdp:
            lines_mdp = cur_mdp.readlines()
            nr_mdp_lines = len(lines_mdp)
            for i in range(nr_mdp_lines):
                spl = lines_mdp[i].split(' ')
                if spl[0] != parameter:
                    with open(new, 'a') as new_mdp:
                        new_mdp.write(lines_mdp[i])
                else:
                    spl[-1] = str(value) + '\n'
                    new_spl = ' '.join(spl)
                    with open(new, 'a') as new_mdp:
                        new_mdp.write(new_spl)
    
    if test == 'no':
        
        mdp_remove_command = 'rm ' + old
        mdp_name_change_command = 'mv ' + new + ' ' + old                   
    
        subprocess.call(mdp_remove_command, shell=True)                   
        subprocess.call(mdp_name_change_command, shell=True)
    
def xtc_maker(all_atoms, lipids, shock_nr, temp_gro_name, rp, test='no'):
    """
    Generates xtc files from gro files each pumping cycle in order to be able
    to concatenate all files at the end of the process.

    Parameters:
    -----------
    lipids: MDAnalysis atomgroup of the lipids in the system
    shock_nr: integer
        shock number

    Returns:
    --------
    generates an xtc file
    """
    if rp == 'no':
        lipids.write(temp_gro_name)
    else:
        all_atoms.write(temp_gro_name)

    if len(str(shock_nr)) == 1:
        new_name = 'vesicle_sA' + str(shock_nr) + '_t.xtc'
    elif len(str(shock_nr)) == 2:
        new_name = 'vesicle_sB' + str(shock_nr) + '_t.xtc'
    else:
        new_name = 'vesicle_sC' + str(shock_nr) + '_t.xtc'

    if test=='no':
        
        convcommand = 'gmx trjconv -f ' + temp_gro_name + ' -o ' + new_name
        subprocess.call(convcommand, shell=True)
    
        rmcommand = 'rm ' + temp_gro_name
        subprocess.call(rmcommand, shell=True)
        
    else:
        return new_name


    
         