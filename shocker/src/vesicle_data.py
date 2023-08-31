#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 18:02:58 2022

@author: marco
"""
import numpy as np
import MDAnalysis.analysis.leaflet
import pyvista as pv
from numpy import linalg as LA


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


def concentration(volume, nr_ions):
    """
    Calculates the ion concentration in mol/L from a given volume and nr of
    ions.

    Parameters:
    -----------
    volume: float
        closed volume the ions reside in
    nr_ions: int
        number of ions in the volume

    Returns:
    --------
    ion concentration
    """
    mol = nr_ions/6.022e23
    lit = volume*1e-24

    concentration = mol/lit

    return concentration


def shape(V0, Ad0, vol, inner, outer):
    """
    Calculates the reduced volume and the reduced area difference

    Parameters
    ----------
    V0 : float
        Initial volume of the vesicle.
    Ad0 : float
        Initial area difference.
    volume : float
        Current vesicle volume.
    inner : float
        current inner vesicle area.
    outer : float
        current outer vesicle area.

    Returns
    -------
    The reduced volume and the reduced area difference.
    """
    rv = vol/V0
    ad = outer - inner
    rad = ad/Ad0

    return rv, rad


def gyration_tensor(cloud):
    """
    Calculates the gyration tensor of a set of points in 3D space

    Parameters
    ----------
    cloud : array
        positions of particle type.

    Returns
    -------
    Gyration tensor
    """
    com = center_of_geometry(cloud)
    corrected_cloud = cloud - com
    N = len(corrected_cloud)

    g_tensor = np.zeros((3, 3))

    for i in corrected_cloud:
        g_tensor[0][0] = g_tensor[0][0] + (i[0]*i[0])
        g_tensor[0][1] = g_tensor[0][1] + (i[0]*i[1])
        g_tensor[0][2] = g_tensor[0][2] + (i[0]*i[2])
        g_tensor[1][0] = g_tensor[1][0] + (i[1]*i[0])
        g_tensor[1][1] = g_tensor[1][1] + (i[1]*i[1])
        g_tensor[1][2] = g_tensor[1][2] + (i[1]*i[2])
        g_tensor[2][0] = g_tensor[2][0] + (i[2]*i[0])
        g_tensor[2][1] = g_tensor[2][1] + (i[2]*i[1])
        g_tensor[2][2] = g_tensor[2][2] + (i[2]*i[2])

    g_tensor = g_tensor/N

    return g_tensor


def delta(tensor):
    """
    Calculates the sphericity tensor

    Parameters
    ----------
    tensor : array
        gyration tensor.

    Returns
    -------
    Sphericity tensor delta.
    """
    values, vectors = LA.eig(tensor)

    trace = np.matrix.trace(tensor)
    lab_bar = trace/3

    c = 0
    for i in values:
        c = c + ((i - lab_bar)**2)

    d = (c/((trace)**2))*(3/2)

    return d


class VolArea():
    """
    Class for calculating the volume and area of a given lipid vesicle
    according to a triangulated surface.

    Parameters:
    -----------
    lipids: MDAnalysis atomgroup
        all lipids in the system
    box_dim: array (3,)
        dimensions of the simulation box
    selection: string
        name of the atom type used to generate a triangulated surface
    """

    def __init__(self,
                 all_atoms,
                 box_dim,
                 selection):

        self.all_atoms = all_atoms
        self.box_dim = box_dim
        self.selection = selection

    def boundary_corrector(self):
        """
        Makes a structure whole again once it passes through the periodic
        boundaries. Only the GL1 beads are targeted since these are used for
        generating a triangulated surface of the vesicle leaflets. The center
        of geometry of the smaller cut-off vesicle parts is calculated and
        according to the distance to the nearest box edge the vesicle part is
        relocated to the correct position in order to construct a complete
        pointcloud.

        Returns:
        --------
        MDAnalysis atomgroup of the GL1 beads in the system, as a whole point
        cloud.
        """
        cloud_particles = MDAnalysis.analysis.leaflet.LeafletFinder(
            self.all_atoms,
            self.selection)

        group_len = []
        for i, _ in enumerate(cloud_particles.groups()):
            group_len.append(len(cloud_particles.group(i)))
        list.sort(group_len)

        base_index = []
        for i, _ in enumerate(cloud_particles.groups()):
            if len(cloud_particles.group(i)) >= group_len[-2]:
                base_index.append(i)

        base = cloud_particles.group(base_index[0]) + \
            cloud_particles.group(base_index[1])

        for i in cloud_particles.groups():
            if len(i) < group_len[-2]:

                cluster = i.positions

                cog = center_of_geometry(cluster)

                absvalues = list(abs(self.box_dim[:3] - cog)) +\
                    list(abs([0, 0, 0] - cog))

                minimum = np.min(absvalues)

                min_index = list(absvalues).index(minimum)

                if min_index < 3:
                    for particle in cluster:
                        particle[min_index] = particle[min_index] -\
                            self.box_dim[min_index]

                if min_index > 2:
                    for particle in cluster:
                        particle[min_index-3] = particle[min_index-3] +\
                            self.box_dim[min_index-3]

                i.positions = cluster

                base = base + i
        return base

    def volume_area(self):
        """
        Uses a triangulated surface to calculate the volume of the vesicle and
        the surface area of both leaflets. If a structure is cut off by
        periodic box boundaries, the structure is made whole again first with
        'boundary corrector'.

        Returns:
        --------
        the area of the inner leaflet, the area of the outer leaflet and the
        vesicle volume.
        """
        leaflets = MDAnalysis.analysis.leaflet.LeafletFinder(self.all_atoms,
                                                             self.selection)
        outer = leaflets.groups(0).positions
        inner = leaflets.groups(1).positions

        # if len(leaflets.groups()) > 2:

        #     pointcloud = self.boundary_corrector()
        #     leaflets = MDAnalysis.analysis.leaflet.LeafletFinder(
        #         pointcloud,
        #         self.selection)

        #     outer = leaflets.groups(0).positions
        #     inner = leaflets.groups(1).positions

        result = None
        while result is None:
            try:
                outer_mesh = mesh_maker(outer)
                outer_area = outer_mesh.area/100
                inner_mesh = mesh_maker(inner)
                inner_area = inner_mesh.area/100

                volume_in = inner_mesh.volume/1000
                box_vol = np.prod(self.box_dim[:3])/1000
                vol_ves_out = outer_mesh.volume/1000
                volume_out = box_vol - vol_ves_out
                result = volume_out
            except RuntimeError:
                pass

        return inner_area, outer_area, volume_in, volume_out


class IonConcentration():
    """
    Class for calculating the ion concentration in the inner and outer
    compartment.
    """

    def __init__(self,
                 bin_size,
                 box_dim,
                 w_clusters):

        self.bin_size = bin_size
        self.box_dim = box_dim
        self.w_clusters = w_clusters
        self.nr_bins = [int(x/self.bin_size) for x in self.box_dim]

    def bin_converter_i(self, target_atoms):
        """
        Each element of a set of particles with coordinates 'target_atoms' in a
        simulationbox of size 'box_dim' is assigned to a bin of size
        'box_dim/nr_bins'. Units are in Angström.

        Parameters:
        -----------
        target_atoms: array
            coordinates of particles in a simulation box

        Returns:
        --------
        array (x,3) containing the positions of the bins containing the
        particles in 'target_atoms'
        """
        x_pos = np.floor(target_atoms[:, 0]/self.bin_size)
        y_pos = np.floor(target_atoms[:, 1]/self.bin_size)
        z_pos = np.floor(target_atoms[:, 2]/self.bin_size)

        binsfloat = np.stack((x_pos, y_pos, z_pos), axis=-1)
        bins = np.int_(binsfloat)

        for i in bins:
            for pos in range(3):
                if i[pos] > self.nr_bins[pos]:
                    i[pos] = i[pos] - self.nr_bins[pos]
                if i[pos] < 0:
                    i[pos] = self.nr_bins[pos] + i[pos]
        return bins

    def ion_counter(self, i_all):
        """
        Finds the indices of the clustered water bins 'w_cluster' in the
        complete list of water bins 'w_all'. The search continues until the
        desired number of indices 'nr_remove' is found. In other words: it
        finds the indices of the water particles residing in the clustered bins

        Parameters:
        -----------
        w_all: array (x, 3)
            a list in which each water particle coordinate in the
            system is converted to the position of the bin it resides in.
        w_cluster: array (x, 3)
            a list of positions of clustered water-bins (vesicle interior)

        Returns:
        --------
        array, indices of the water particles we want to remove from the
        vesicle interior according to the total water list
        """
        nr_ions = []

        for cluster in self.w_clusters:

            indices = []
            for i, _ in enumerate(cluster):

                bin_index = np.where((i_all[:, 0] == cluster[i][0]) &
                                     (i_all[:, 1] == cluster[i][1]) &
                                     (i_all[:, 2] == cluster[i][2]))[0]

                for index in bin_index:
                    indices.append(int(index))
            nr_ions.append(len(indices))

        return nr_ions
