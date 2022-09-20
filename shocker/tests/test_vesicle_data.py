# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 11:41:50 2022

@author: s103567
"""

import numpy as np
import os
import MDAnalysis.analysis.leaflet
import MDAnalysis as mda
import pytest
import shocker.src.vesicle_data as vd
from shocker.src.vesicle_data import VolArea
from shocker.src.vesicle_data import IonConcentration
from shocker import TEST_DATA

data_path = TEST_DATA


def test_center_of_geometry():

    pointcloud = [(1, 2, 3), (4, 3, 2), (3, 2, 4)]
    cog_result = np.array([2.66666667, 2.33333333, 3.])
    assert all(vd.center_of_geometry(pointcloud)) == all(cog_result)


def test_concentration():

    volume = 13000
    nr_ions = 10000
    result_conc = 1.28
    function_result_conc = round(vd.concentration(volume, nr_ions), 2)
    assert result_conc == function_result_conc


def test_shape():

    V0 = 13500
    Ad0 = 1200
    volume = 13000
    inner = 2780
    outer = 3925

    result_rv = 0.96
    result_rad = 0.95

    function_result = vd.shape(V0, Ad0, volume, inner, outer)
    function_result_rv = round(function_result[0], 2)
    function_result_rad = round(function_result[1], 2)
    assert result_rv == function_result_rv
    assert result_rad == function_result_rad


def test_gyration_tensor():

    pointcloud = [(4, 4, 4), (3, 2, 3), (2, 4, 3)]
    result = np.array([[0.66666667, 0.        , 0.33333333],
                       [0.        , 0.88888889, 0.22222222],
                       [0.33333333, 0.22222222, 0.22222222]])
    function_result = vd.gyration_tensor(pointcloud)
    assert result.all() == function_result.all()


def test_delta():

    tensor = np.array([[0.66666667, 0.        , 0.33333333],
                       [0.        , 0.88888889, 0.22222222],
                       [0.33333333, 0.22222222, 0.22222222]])

    result = 0.26
    function_result = round(vd.delta(tensor), 2)
    assert result == function_result


@pytest.fixture
def areavolume_object():

    file_path = os.path.join(data_path, 'vesicle_data_test.gro')
    universe = mda.Universe(file_path)
    all_atoms = universe.select_atoms('not resname W')
    tj = universe.trajectory[0]
    box_dim = tj.dimensions
    selection = 'name GL1'

    kwargs = {"all_atoms": all_atoms,
              "box_dim": box_dim,
              "selection": selection}

    return VolArea(**kwargs)


def test_volume_area(areavolume_object):

    result_area_in = 2282
    result_area_out = 3597
    result_volume = 8673

    vesicle_data = areavolume_object.volume_area()

    function_result_area_in = vesicle_data[0]
    function_result_area_out = vesicle_data[1]
    function_result_volume = vesicle_data[2]

    assert result_area_in == int(function_result_area_in)
    assert result_area_out == int(function_result_area_out)
    assert result_volume == int(function_result_volume)


@pytest.fixture
def ionconcentration_object():

    file_path = os.path.join(data_path, 'vesicle_data_test.gro')
    universe = mda.Universe(file_path)
    tj = universe.trajectory[0]
    box_dim = tj.dimensions
    bin_size = 10
    w_clusters = np.array([[[1, 2, 3], [2, 3, 4],
                            [5, 4, 3], [5, 5, 5],
                            [4, 4, 4]]])

    kwargs = {"bin_size": bin_size,
              "box_dim": box_dim,
              "w_clusters": w_clusters}

    return IonConcentration(**kwargs)


def test_bin_converter_i(ionconcentration_object):

    test_system = np.array([[22, 44, 55], [102, 35, 250],
                            [75, 32, 11], [90, 23, 30]])
    result = np.array([[2, 4, 5], [10, 3, 25],
                       [7, 3, 1], [9, 2, 3]])
    function_result = ionconcentration_object.bin_converter_i(test_system)

    for i, _ in enumerate(result):
        assert all(function_result[i] == result[i])


def test_ion_counter(ionconcentration_object):

    i_all = np.array([[1, 2, 3], [2, 3, 4]])
    result = [2]
    function_result = ionconcentration_object.ion_counter(i_all)

    assert result == function_result
