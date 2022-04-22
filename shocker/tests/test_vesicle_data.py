# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 11:41:50 2022

@author: s103567
"""

import numpy as np
import MDAnalysis.analysis.leaflet
import MDAnalysis as mda
import pytest
import shocker.src.vesicle_data as vd
from shocker.src.vesicle_data import VolArea



def test_center_of_geometry():
    
    pointcloud = [(1,2,3),(4,3,2),(3,2,4)]
    cog_result = np.array([2.66666667, 2.33333333, 3.        ])
    assert all(vd.center_of_geometry(pointcloud)) == all(cog_result)
    
@pytest.fixture
def areavolume_object():
    
    universe = mda.Universe('test_data/vesicle_data_test.gro')                
    lipids = universe.select_atoms('not resname W')
    tj = universe.trajectory[0]
    box_dim = tj.dimensions
    selection = 'name GL1'
    
    kwargs = {"lipids": lipids, "box_dim": box_dim, "selection": selection}

    return VolArea(**kwargs)

def test_volume_area(areavolume_object):
    
    result_area_in = 228235
    result_area_out = 359754
    result_volume = 8673499
    
    vesicle_data = areavolume_object.volume_area()
    
    function_result_area_in = vesicle_data[0]
    function_result_area_out = vesicle_data[1]
    function_result_volume = vesicle_data[2]
    
    assert result_area_in == int(function_result_area_in)
    assert result_area_out == int(function_result_area_out)
    assert result_volume == int(function_result_volume)
    
    
    
    