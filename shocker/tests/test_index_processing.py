#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 09:19:35 2022

@author: marco
"""

import os
import shocker.src.index_processing as ip
import MDAnalysis as mda
from shocker import TEST_DATA

data_path = TEST_DATA


def test_ag_extractor():

    file_path = os.path.join(data_path, 'index_test.gro')
    universe = mda.Universe(file_path)
    all_atoms = universe.select_atoms('all')
    index_file = os.path.join(data_path, 'index.ndx')

    result = ([universe.select_atoms('resname POPC'),
               universe.select_atoms('resname W')],
              ['membrane', 'water'])

    function_result = ip.ag_extractor(all_atoms, index_file)

    assert result == function_result


def test_ag_finder():

    file_path = os.path.join(data_path, 'index_test.gro')
    universe = mda.Universe(file_path)
    all_atoms = universe.select_atoms('all')
    index_file = os.path.join(data_path, 'index.ndx')
    groups_and_names = ip.ag_extractor(all_atoms, index_file)

    result = universe.select_atoms('resname W')

    function_result = ip.ag_finder('water', groups_and_names)

    assert result == function_result


def test_name_concatenator():

    name_list = ['membrane', 'water']

    result = 'membrane water'

    function_result = ip.name_concatenator(name_list)

    assert result == function_result
