# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 14:04:27 2022

@author: s103567
"""
import os
import shocker.src.file_mod as fm
import numpy as np
import subprocess
import MDAnalysis as mda
from shocker import TEST_DATA

data_path = TEST_DATA


def test_name_generator():

    old_name = 'testname.gro'
    shock_nr = 44

    result = 'testname_s44.gro'
    result_function = fm.name_generator(old_name, shock_nr, 'yes')

    assert result_function == result


def test_mdp_value_changer():

    mdp_old = os.path.join(data_path, 'test.mdp')

    parameter = 'nsteps'
    value = 50000

    result = os.path.join(data_path, 'result.mdp')
    function_result = os.path.join(data_path, 'new.mdp')

    fm.mdp_value_changer(mdp_old, function_result, parameter, value, 'yes')

    with open(result, 'r') as res:
        reslines = res.readlines()

    with open(function_result, 'r') as funcres:
        funcreslines = funcres.readlines()

    assert reslines == funcreslines

    rm_command = 'rm {}'.format(function_result)
    subprocess.call(rm_command, shell=True)


def test_xtc_maker():

    gro_file = os.path.join(data_path, 'water_remover_test.gro')
    universe = mda.Universe(gro_file)
    lipids = universe.select_atoms('not resname W')
    shock_nr = 44

    temp_gro_name = os.path.join(data_path, 'only_lipids_temp.gro')
    temp_gro_name_result = os.path.join(data_path, 'only_lipids_result.gro')

    new_xtc_name_result = 'vesicle_sB44_t.xtc'
    new_xtc_name_function_result = fm.xtc_maker(lipids, shock_nr,
                                                temp_gro_name, 'yes')

    assert new_xtc_name_result == new_xtc_name_function_result

    with open(temp_gro_name, 'r') as res:
        reslines = res.readlines()

    with open(temp_gro_name_result, 'r') as funcres:
        funcreslines = funcres.readlines()

    assert reslines == funcreslines


def test_gro_to_np():

    gro_file = os.path.join(data_path, 'gro_to_np_test.gro')
    result = np.array([[39.29, 25.47, 15.02],
                       [21.05, 28.94,  3.17],
                       [13.11, 38.59, 21.29],
                       [32.91, 36.68,  3.39],
                       [4.00, 4.00, 4.00],
                       [3.17, 22.86, 0.42],
                       [1.00, 1.00, 1.00]])

    function_result = fm.gro_to_np(gro_file)

    for i, v in enumerate(result):
        for j, _ in enumerate(v):
            assert round(function_result[i][j], 2) == result[i][j]
