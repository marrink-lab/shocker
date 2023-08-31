# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 13:44:26 2022

@author: s103567
"""
import pytest
import numpy as np
import shocker.src.water_repositioning as wr
import MDAnalysis as mda
import os
from shocker.src.water_repositioning import Mover
from shocker import TEST_DATA

data_path = TEST_DATA


def bin_converter(target_atoms):

    x_pos = np.floor(target_atoms[:, 0]/10)
    y_pos = np.floor(target_atoms[:, 1]/10)
    z_pos = np.floor(target_atoms[:, 2]/10)

    binsfloat = np.stack((x_pos, y_pos, z_pos), axis=-1)
    bins = np.int_(binsfloat)

    return bins


def test_index_finder():

    test_data_all = np.array([[2, 3, 4], [2, 2, 1], [6, 4, 3],
                              [1, 1, 1], [2, 2, 2], [2, 1, 2],
                              [4, 3, 2], [8, 7, 7], [9, 8, 7]])
    test_data_bin = np.array([2, 2, 2])

    result = [4]
    function_result = wr.index_finder(test_data_all, test_data_bin)

    assert result == function_result


@pytest.fixture
def mover_object():

    receiving_cluster = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0],
                                  [0, 1, 0], [1, 1, 0], [2, 1, 0],
                                  [0, 2, 0], [1, 2, 0], [2, 2, 0],
                                  [0, 0, 1], [1, 0, 1], [2, 0, 1],
                                  [0, 1, 1], [1, 1, 1], [2, 1, 1],
                                  [0, 2, 1], [1, 2, 1], [2, 2, 1],
                                  [0, 0, 2], [1, 0, 2], [2, 0, 2],
                                  [0, 1, 2], [1, 1, 2], [2, 1, 2],
                                  [0, 2, 2], [1, 2, 2], [2, 2, 2]])

    file_path = os.path.join(data_path, 'test_data_water_move.gro')
    universe = mda.Universe(file_path)
    all_atoms = universe.select_atoms('all')
    bin_size = 10
    box_dim = [40, 40, 40]
    nr_removed = 5
    all_bins = bin_converter(all_atoms.positions)
    search_volume = 3

    water_group = universe.select_atoms('resname W')

    kwargs = {"receiving_cluster": receiving_cluster,
              "all_atoms": all_atoms,
              "bin_size": bin_size,
              "nr_removed": nr_removed,
              "box_dim": box_dim,
              "all_bins": all_bins,
              "water_group": water_group,
              "search_volume": search_volume}

    return Mover(**kwargs)


@pytest.fixture
def mover_object_aa():

    receiving_cluster = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0],
                                  [0, 1, 0], [1, 1, 0], [2, 1, 0],
                                  [0, 2, 0], [1, 2, 0], [2, 2, 0],
                                  [0, 0, 1], [1, 0, 1], [2, 0, 1],
                                  [0, 1, 1], [1, 1, 1], [2, 1, 1],
                                  [0, 2, 1], [1, 2, 1], [2, 2, 1],
                                  [0, 0, 2], [1, 0, 2], [2, 0, 2],
                                  [0, 1, 2], [1, 1, 2], [2, 1, 2],
                                  [0, 2, 2], [1, 2, 2], [2, 2, 2]])

    file_path = os.path.join(data_path, 'test_data_water_move_aa.gro')
    universe = mda.Universe(file_path)
    all_atoms = universe.select_atoms('all')

    bin_size = 10
    box_dim = [10, 10, 23]
    nr_removed = 5
    all_bins = bin_converter(all_atoms.positions)
    search_volume = 3

    water_group = universe.select_atoms('resname TIP3')

    kwargs = {"receiving_cluster": receiving_cluster,
              "all_atoms": all_atoms,
              "bin_size": bin_size,
              "nr_removed": nr_removed,
              "box_dim": box_dim,
              "all_bins": all_bins,
              "water_group": water_group,
              "search_volume": search_volume}

    return Mover(**kwargs)


def test_repositioning_bin_identifier(mover_object):

    resulting_bins = mover_object.repositioning_bin_identifier()

    assert len(resulting_bins) == mover_object.nr_removed

    doubles = 0
    bin_list = [list(x) for x in resulting_bins]
    for i in bin_list:
        if bin_list.count(i) > 1:
            doubles = doubles + 1
    assert doubles == 0


def test_position_generator(mover_object):

    chosen_bins = np.array([[3, 0, 0],
                            [3, 3, 1],
                            [2, 1, 3],
                            [0, 0, 3],
                            [3, 3, 3]])

    result_min_dist = 3
    function_result = mover_object.position_generator(chosen_bins)

    assert function_result[1] > result_min_dist

    assert len(function_result[0]) == mover_object.nr_removed


def test_position_generator_aa(mover_object_aa):

    chosen_bins = np.array([[6, 4, 12],
                            [1, 3, 8],
                            [4, 9, 14],
                            [7, 0, 8],
                            [6, 5, 8]])
    oxy_pos = 0

    result_min_dist = 3
    function_result = mover_object_aa.position_generator_aa(chosen_bins,
                                                            oxy_pos)

    assert function_result[1] > result_min_dist

    assert len(function_result[0]) == mover_object_aa.nr_removed


def test_water_repositioning_gro(mover_object):

    indices = [4, 6, 8]
    new_pos = np.array([[4, 4, 4],
                        [1, 1, 1],
                        [2, 2, 2]])

    result = os.path.join(data_path, 'result_data_water_move.gro')
    function_result = os.path.join(data_path, 'function_result_water_move.gro')

    mover_object.water_repositioning_gro(indices, new_pos, function_result)

    with open(result, 'r') as res:
        reslines = res.readlines()

    with open(function_result, 'r') as funcres:
        funcreslines = funcres.readlines()

    assert reslines == funcreslines


def test_water_repositioning_gro_aa(mover_object_aa):

    indices = [3, 6, 9]
    new_pos = np.array([[4, 4, 4],
                        [1, 1, 1],
                        [2, 2, 2]])

    result = os.path.join(data_path, 'result_data_water_move_aa.gro')
    function_result = os.path.join(data_path,
                                   'function_result_water_move_aa.gro')

    mover_object_aa.water_repositioning_gro_aa(indices,
                                               new_pos,
                                               function_result,
                                               0)

    with open(result, 'r') as res:
        reslines = res.readlines()

    with open(function_result, 'r') as funcres:
        funcreslines = funcres.readlines()

    assert reslines == funcreslines
