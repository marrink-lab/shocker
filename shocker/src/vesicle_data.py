#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 18:02:58 2022

@author: marco
"""
import subprocess
import numpy as np
import MDAnalysis.analysis.leaflet
import pyvista as pv

def mesh_maker(cloud):
    """
    Generates a triangulated mesh of a pointcloud

    Parameters:
    -----------
    cloud: array (x,3)
        array of coordinates in space

    Returns:
    --------
    :class: core.pointset.PolyData
        a triangulated surface
    """
    points = pv.wrap(cloud)
    surf = points.reconstruct_surface(nbr_sz=6)

    return surf

def center_of_geometry(points):
    """
    Calculates the center of geometry of a set of coordinates and stores it
    in a numpy array.

    Parameters:
    -----------
    points: array (x,3)
        a set of coordinates

    Returns:
    --------
    array (3,) containing the center of geometry of 'points'
    """
    nr_points = len(points)
    x_sum = 0
    y_sum = 0
    z_sum = 0

    for i in points:
        x_sum = x_sum + i[0]
        y_sum = y_sum + i[1]
        z_sum = z_sum + i[2]

    cog = np.array([x_sum/nr_points, y_sum/nr_points, z_sum/nr_points])

    return cog

class VolArea():
    """
    Class for calculating the volume and area of a given lipid vesicle according
    to a triangulated surface
    
    Parameters:
    -----------
    lipids: MDAnalysis atomgroup
        all lipids in the system
    box_dim: array (3,)
        dimensions of the simulation box
    selection: string
        name of the atom type used to generate a triangulated surface
    """
    
    def __init__(self, lipids, box_dim, selection):
        
        self.lipids = lipids
        self.box_dim = box_dim
        self.selection = selection
        
    def boundary_corrector(self):
        """
        Makes a structure whole again once it passes through the periodic boundaries.
        Only the GL1 beads are targeted since these are used for generating a
        triangulated surface of the vesicle leaflets. The center of geometry of the
        smaller cut-off vesicle parts is calculated and according to the distance to the
        nearest box edge the vesicle part is relocated to the correct position in order
        to construct a complete pointcloud.
    
        Returns:
        --------
        MDAnalysis atomgroup of the GL1 beads in the system, as a whole point cloud
        """
        cloud_particles = MDAnalysis.analysis.leaflet.LeafletFinder(self.lipids, self.selection)
    
        group_len = []
        for i, _ in enumerate(cloud_particles.groups()):
            group_len.append(len(cloud_particles.group(i)))
        list.sort(group_len)
    
        base_index = []
        for i, _ in enumerate(cloud_particles.groups()):
            if len(cloud_particles.group(i)) >= group_len[-2]:
                base_index.append(i)
        base = cloud_particles.group(base_index[0]) + cloud_particles.group(base_index[1])
    
        for i in cloud_particles.groups():
            if len(i) < group_len[-2]:
    
                cluster = i.positions
    
                cog = center_of_geometry(cluster)
    
                absvalues = list(abs(self.box_dim[:3] - cog)) + list(abs([0,0,0] - cog))
    
                minimum = np.min(absvalues)
    
                min_index = list(absvalues).index(minimum)
    
                if min_index < 3:
                    for particle in cluster:
                        particle[min_index] = particle[min_index] - self.box_dim[min_index]
    
                if min_index > 2:
                    for particle in cluster:
                        particle[min_index-3] = particle[min_index-3] + self.box_dim[min_index-3]
    
                i.positions = cluster
    
                base = base + i
        return base
    
    def volume_area(self):
        """
        Uses a triangulated surface to calculate the volume of the vesicle and the
        surface area of both leaflets. If a structure is cut off by periodic box
        boundaries, the structure is made whole again first with 'boundary corrector'.
        
        Returns:
        --------
        the area of the inner leaflet, the area of the outer leaflet and the vesicle
        volume.
        """
        leaflets = MDAnalysis.analysis.leaflet.LeafletFinder(self.lipids, self.selection)
        outer = leaflets.groups(0).positions
        inner = leaflets.groups(1).positions
    
        if len(leaflets.groups()) > 2:
    
            pointcloud = self.boundary_corrector()
            leaflets = MDAnalysis.analysis.leaflet.LeafletFinder(pointcloud, self.selection)
            outer = leaflets.groups(0).positions
            inner = leaflets.groups(1).positions
    
        result = None
        while result is None:
            try:
                outer_mesh = mesh_maker(outer)
                outer_area = outer_mesh.area
                inner_mesh = mesh_maker(inner)
                inner_area = inner_mesh.area
    
                volume = inner_mesh.volume
                result = volume
            except RuntimeError:
                pass
    
        return inner_area, outer_area, volume


