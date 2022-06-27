# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 11:39:21 2022

@author: s103567
"""

import numpy as np
import os
import pytest
import shocker.src.cluster as cl
import MDAnalysis as mda
from shocker.src.cluster import Cluster
from shocker import TEST_DATA

data_path = TEST_DATA


@pytest.fixture
def multiple_bin_list():
    mbl = np.array([[1, 2, 3], [4, 4, 4], [5, 4, 3],
                    [4, 4, 4], [7, 4, 3]])
    return mbl


def test_multiples_remover(multiple_bin_list):
    bins = multiple_bin_list
    result = np.array([[1, 2, 3], [4, 4, 4],
                       [5, 4, 3], [7, 4, 3]])
    function_result = cl.multiples_remover(bins)

    for i, _ in enumerate(result):
        assert all(function_result[i] == result[i])


def test_zero_finder():
    test_system = np.ones((4, 4, 4))
    test_system[2][2][2] = 0
    result = (2, 2, 2)
    function_result = cl.zero_finder(test_system)

    assert function_result == result


@pytest.fixture
def cluster_object():
    file_path = os.path.join(data_path, 'vesicle_system.gro')
    universe = mda.Universe(file_path)
    all_atoms = universe.select_atoms('all')
    tj = universe.trajectory[0]
    box_dim = tj.dimensions
    bin_size = 10

    kwargs = {"all_atoms": all_atoms, "box_dim": box_dim, "bin_size": bin_size}

    return Cluster(**kwargs)


@pytest.fixture
def cluster_object_aa():
    file_path = os.path.join(data_path, 'doublebilayer_test.gro')
    universe = mda.Universe(file_path)
    all_atoms = universe.select_atoms('all')
    tj = universe.trajectory[0]
    box_dim = tj.dimensions
    bin_size = 10

    kwargs = {"all_atoms": all_atoms,
              "box_dim": box_dim,
              "bin_size": bin_size}

    return Cluster(**kwargs)


def test_bin_converter_l(cluster_object):
    test_system = np.array([[22, 44, 55], [102, 35, 250],
                            [75, 32, 11], [90, 23, 30]])
    result = np.array([[2, 4, 5], [10, 3, 25],
                       [7, 3, 1], [9, 2, 3]])
    function_result = cluster_object.bin_converter_l(test_system)

    for i, _ in enumerate(result):
        assert all(function_result[i] == result[i])


def test_bin_system_maker(cluster_object):
    test_system = np.array([[2, 2, 2], [0, 2, 1],
                            [1, 2, 2], [0, 0, 0]])
    result = np.zeros((cluster_object.nr_bins[0],
                       cluster_object.nr_bins[1],
                       cluster_object.nr_bins[2]))
    result[2][2][2] = 1
    result[0][2][1] = 1
    result[1][2][2] = 1
    result[0][0][0] = 1

    function_result = cluster_object.bin_system_maker(test_system)

    for i, v in enumerate(result):
        for j, _ in enumerate(v):
            assert all(function_result[i][j] == result[i][j])


def test_neighbor_view(cluster_object):
    function_storage = []
    cur = [1, 1, 1]
    direction = (1, 0, 0)
    function_system = np.zeros((3, 3, 3))
    result_storage = [(2, 1, 1)]
    result_system = np.array([[[0., 0., 0.],
                               [0., 0., 0.],
                               [0., 0., 0.]],

                              [[0., 0., 0.],
                               [0., 0., 0.],
                               [0., 0., 0.]],

                              [[0., 0., 0.],
                               [0., -1., 0.],
                               [0., 0., 0.]]])
    cluster_object.neighbor_view(direction,
                                 function_system,
                                 function_storage,
                                 cur)

    assert function_storage == result_storage
    for i, v in enumerate(result_system):
        for j, _ in enumerate(v):
            assert all(function_system[i][j] == result_system[i][j])


def test_cluster_finder(cluster_object):

    out_cluster_test = []
    for x in range(cluster_object.nr_bins[0]):
        for y in range(cluster_object.nr_bins[1]):
            for z in range(cluster_object.nr_bins[2]):
                out_cluster_test.append((x, y, z))

    test_system = np.zeros((cluster_object.nr_bins[0],
                            cluster_object.nr_bins[1],
                            cluster_object.nr_bins[2]))
    for x in range(20, 23):
        for y in range(20, 23):
            for z in range(20, 23):
                test_system[x][y][z] = 1
                out_cluster_test.remove((x, y, z))
    test_system[21][21][21] = 0

    in_cluster_test = [(21, 21, 21)]

    cluster_function = cluster_object.cluster_finder(test_system)

    assert set(cluster_function[0]) == set(out_cluster_test)
    assert set(cluster_function[1]) == set(in_cluster_test)
