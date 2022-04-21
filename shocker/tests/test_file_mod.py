# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 14:04:27 2022

@author: s103567
"""
import pytest
import shocker.src.file_mod as fm
import subprocess
import MDAnalysis as mda

def test_name_generator():
    
    old_name = 'testname.gro'
    shock_nr = 44
    
    result = 'testname_s44.gro'
    result_function = fm.name_generator(old_name, shock_nr, 'yes')
    
    assert result_function == result
    
def test_mdp_value_changer():
    
    mdp_old = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/test.mdp'
    parameter = 'nsteps'
    value = 50000
    
    result = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/result.mdp'
    function_result = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/new.mdp'
    
    fm.mdp_value_changer(mdp_old, function_result, parameter, value, 'yes')
    
    with open(result, 'r') as res:
        reslines = res.readlines()
    
    with open(function_result, 'r') as funcres:
        funcreslines = funcres.readlines()
        
    assert reslines == funcreslines
    
    rm_command = 'rm {}'.format(function_result)
    subprocess.call(rm_command, shell=True)
    
def test_xtc_maker():
    
    gro_file = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/water_remover_test.gro'
    universe = mda.Universe(gro_file)
    lipids = universe.select_atoms('not resname W')
    shock_nr = 44
    
    temp_gro_name = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/only_lipids_temp.gro'
    temp_gro_name_result = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/only_lipids_result.gro'
    
    new_xtc_name_result = 'vesicle_sB44_t.xtc'
    new_xtc_name_function_result = fm.xtc_maker(lipids, shock_nr, temp_gro_name, 'yes')
    
    assert new_xtc_name_result == new_xtc_name_function_result
    
    with open(temp_gro_name, 'r') as res:
        reslines = res.readlines()
    
    with open(temp_gro_name_result, 'r') as funcres:
        funcreslines = funcres.readlines()
        
    assert reslines == funcreslines
    
    