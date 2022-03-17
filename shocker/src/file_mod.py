#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 18:02:38 2022

@author: marco
"""
import subprocess

def mdp_value_changer(old, parameter, value):
    """
    Modifies a value in the mdp file by writing a new file and deleting the old
    one.
    
    Parameters:
    -----------
    old: string
        path to the old mdp file
    parameter: string
        parameter name to be modified
    value: integer
        new parameter value
        
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
                    with open('new.mdp', 'a') as new_mdp:
                        new_mdp.write(lines_mdp[i])
                else:
                    spl[-1] = str(value) + '\n'
                    new_spl = ' '.join(spl)
                    with open('new.mdp', 'a') as new_mdp:
                        new_mdp.write(new_spl)

    mdp_remove_command = 'rm ' + old
    mdp_name_change_command = 'mv ' + 'new.mdp ' + old                   

    subprocess.call(mdp_remove_command, shell=True)                   
    subprocess.call(mdp_name_change_command, shell=True)
    
    