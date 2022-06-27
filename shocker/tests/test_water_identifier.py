# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 14:18:34 2022

@author: s103567
"""

import numpy as np
import shocker.src.water_identifier as wi
import pytest
import os
import MDAnalysis as mda
from shocker.src.water_identifier import Identifier
from shocker import TEST_DATA

data_path = TEST_DATA


def test_local_to_global():

    no_water = 500
    local_index = [2, 4, 5, 6, 7, 34, 44, 55]

    result = [502, 504, 505, 506, 507, 534, 544, 555]
    function_result = wi.local_to_global(no_water, local_index)

    assert result == function_result


@pytest.fixture
def identifier_object():

    bin_clusters = [[(0, 0, 0), (2, 2, 2), (43, 43, 44)],
                    [(4, 2, 3), (5, 4, 3), (1, 5, 3)],
                    [(13, 34, 23), (12, 41, 16), (33, 22, 11)]]
    file_path = os.path.join(data_path, 'vesicle_system.gro')
    universe = mda.Universe(file_path)
    tj = universe.trajectory[0]
    box_dim = tj.dimensions
    bin_size = 10
    nr_remove = 10
    water = universe.select_atoms('resname W')

    kwargs = {"bin_clusters": bin_clusters,
              "box_dim": box_dim,
              "bin_size": bin_size,
              "nr_remove": nr_remove,
              "water": water}

    return Identifier(**kwargs)


def test_cluster_selecter(identifier_object):

    outer_cluster_test = [(0, 0, 0), (2, 2, 2), (43, 43, 44)]
    inner_cluster_test = [(4, 2, 3), (5, 4, 3),
                          (1, 5, 3), (13, 34, 23),
                          (12, 41, 16), (33, 22, 11)]

    function_clusters = identifier_object.cluster_selecter()

    assert set(function_clusters[0]) == set(inner_cluster_test)
    assert set(function_clusters[1]) == set(outer_cluster_test)


def test_bin_converter_w(identifier_object):

    test_system = np.array([[22, 44, 55], [102, 35, 250],
                            [75, 32, 11], [90, 23, 30]])

    result = np.array([[2, 4, 5], [10, 3, 25],
                       [7, 3, 1], [9, 2, 3]])
    function_result = identifier_object.bin_converter_w(test_system)

    for i, _ in enumerate(result):
        assert all(function_result[i] == result[i])


def test_index_finder_m(identifier_object):

    cluster_bins = np.array([[1, 2, 3], [3, 2, 1],
                             [5, 4, 5], [4, 4, 4],
                             [7, 3, 4], [6, 5, 4]])

    all_bins = np.array([[1, 2, 3], [1, 2, 3], [4, 4, 4], [5, 4, 5],
                         [4, 4, 4], [5, 4, 5], [1, 1, 1], [4, 4, 4],
                         [6, 5, 4], [7, 3, 4], [7, 3, 4], [2, 2, 2],
                         [6, 6, 6], [7, 3, 4], [1, 2, 3], [3, 2, 1]])

    result = [100608, 100609, 100622,
              100611, 100613, 100610,
              100612, 100615, 100617,
              100618, 100621, 100623,
              100616]

    function_result = identifier_object.index_finder_m(all_bins, cluster_bins)

    assert len(function_result) == identifier_object.nr_remove
    for i in function_result:
        assert i in result


def test_index_finder_m_aa(identifier_object):

    cluster_bins = np.array([[1, 2, 3], [3, 2, 1],
                             [5, 4, 5], [4, 4, 4],
                             [7, 3, 4], [6, 5, 4]])

    all_bins = np.array([[1, 2, 3], [1, 2, 3], [4, 4, 4], [5, 4, 5],
                         [4, 4, 4], [5, 4, 5], [1, 1, 1], [4, 4, 4],
                         [6, 5, 4], [7, 3, 4], [7, 3, 4], [2, 2, 2],
                         [6, 6, 6], [7, 3, 4], [1, 2, 3], [3, 2, 1]])
    oxy_pos = 0
    result = [100608, 100609, 100622,
              100611, 100613, 100610,
              100612, 100615, 100617,
              100618, 100621, 100623,
              100616]

    function_result = identifier_object.index_finder_m_aa(all_bins,
                                                          cluster_bins,
                                                          oxy_pos)

    assert len(function_result) == identifier_object.nr_remove
    for i in function_result:
        assert i in result
