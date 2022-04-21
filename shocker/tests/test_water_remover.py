# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 14:42:42 2022

@author: s103567
"""

import subprocess
import pytest
import MDAnalysis as mda
from shocker.src.water_remover import Remover

@pytest.fixture
def remover_object():
    
    universe = mda.Universe('/mnt/c/TempSim4/shocker/shocker/tests/test_data/water_remover_test.gro')                
    all_atoms = universe.select_atoms('all')
    top_old = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/topol_test.top'
    nr_removed = 5
    gro_file = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/water_remover_function_result.gro'
            
    kwargs = {"top_old": top_old, "all_atoms": all_atoms, "nr_removed": nr_removed, "gro_file": gro_file}

    return Remover(**kwargs)

def test_water_remover_gro(remover_object):
    
    indices = [27, 28, 29, 30, 31]
    result = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/water_remover_result.gro'
    
    remover_object.water_remover_gro(indices)
    function_result = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/water_remover_function_result.gro'
    
    with open(result, 'r') as res:
        reslines = res.readlines()
    
    with open(function_result, 'r') as funcres:
        funcreslines = funcres.readlines()
        
    assert reslines == funcreslines
    
def test_water_remover_top(remover_object):
        
    result = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/topol_result.top'
    
    function_result = '/mnt/c/TempSim4/shocker/shocker/tests/test_data/topol_function_result.top'
    remover_object.water_remover_top(function_result)
        
    with open(result, 'r') as res:
        reslines = res.readlines()
    
    with open(function_result, 'r') as funcres:
        funcreslines = funcres.readlines()
        
    assert reslines == funcreslines
    
    rm_command = 'rm {}'.format(function_result)
    subprocess.call(rm_command, shell=True)
    
    