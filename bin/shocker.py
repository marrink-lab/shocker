#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:47:21 2022

@author: marco
"""

import MDAnalysis as mda
import subprocess
import argparse
import MDAnalysis.analysis.leaflet
import numpy as np
import os
import Shocker.shocker.src.cluster as cl
import Shocker.shocker.src.file_mod as fm
import Shocker.shocker.src.water_replacement as wrep
import Shocker.shocker.src.water_identifier as wi
from Shocker.shocker.src.cluster import Cluster
from Shocker.shocker.src.water_identifier import Identifier
from Shocker.shocker.src.water_remover import Remover
from Shocker.shocker.src.water_replacement import Mover
from Shocker.shocker.src.vesicle_data import VolArea

def main():
    
    parser = argparse.ArgumentParser(description='Osmotic shock simulator', add_help=True)
    parser.add_argument('-in','--input', help='mdp file with input data')
    parser.add_argument('-f','--name', help='name of the files')
    parser.add_argument('-t','--topology', help='topology file')
    parser.add_argument('-r','--removed', type=int, help='nr of removed particles per shock')
    parser.add_argument('-e','--end', type=int, help='nr of the last cycle')
    parser.add_argument('-w','--water', default='resname W', help='name of water particles')
    parser.add_argument('-b','--bins', type=int, default=1, help='bin size for clustering (standard size is 1nm^3)')
    parser.add_argument('-s','--start', type=int, default=0, help='define the starting cycle')
    parser.add_argument('-x','--xtc', type=int, default=0, help='frequency of storing position data (default = 0)')
    parser.add_argument('-lip','--lipids', nargs='+', help='the names of the lipids the membrane consists of')
    parser.add_argument('-it', '--iterations', type=int, default=100000, help='the number of iterations per cycle')
    parser.add_argument('-min','--enermin', default='no', help='energy minimization after each pumping cycle (yes or no)')
    parser.add_argument('-mdpmin', '--emruns', nargs='+', help='the .mdp files used for energy minimization')
    parser.add_argument('-emnr', '--enernr', type=int, default=1, help='number of energy minimization rounds per cycle')
    parser.add_argument('-vd', '--vesdata', default='no', help='whether or not leaflet area and vesicle volume is calculated and stored each cycle')
    parser.add_argument('-cp', '--cloudpoints', default='GL1', help='defines which atoms are used to construct a triangulated surface of the vesicle')
    parser.add_argument('-rp', '--replace', default='no', help='whether or not the waterparticles are replaced or removed')
    parser.add_argument('-shock', '--shocktype', default='hypertonic', help='type of osmotick shock (hypertonic or hypotonic)')
    args = parser.parse_args()
    
    # changing the number of iterations per cycle and frequency of position data 
    # storage (default is 100.000 iterations) if needed
    fm.mdp_value_changer(args.input, 'nsteps', args.iterations)
    fm.mdp_value_changer(args.input, 'nstxout-compressed', args.xtc)
    
    # creating a destination folder to store all generated files
    subprocess.call('mkdir shockfiles', shell=True)
    
    # =========================================================================
    # Defining values and executing simulation run
    # =========================================================================
    
    for shock in range(args.start, args.end):
              
        gro_file = args.name + '.gro'
        tpr_file = args.name + '.tpr'
        gro_file_2 = args.name + '2.gro'
        xtc_file = args.name + '.xtc'
        topology_file = args.topology + '.top'
        topology_file_2 = args.topology + '2.top'
        water_name = args.water
        
        # Energy minimization if needed
        if args.enermin == 'yes':
            for i in range(args.enernr):
                for mdp in args.emruns:
                    
                    grompp_command_em = 'gmx grompp -f ' + mdp + '.mdp' + ' -c ' + gro_file + \
                        ' -p ' + topology_file + ' -maxwarn 1'
                    mdrun_command_em = 'gmx mdrun -v'
                    em_name_change_command = 'mv confout.gro' + ' ' + gro_file
                    tpr_remove_command = 'rm topol.tpr'
                    
                    subprocess.call(grompp_command_em, shell=True)
                    subprocess.call(mdrun_command_em, shell=True)
                    subprocess.call(em_name_change_command, shell=True)
                    subprocess.call(tpr_remove_command, shell=True)
                    
        # Executing the actual simulation 
        restart_condition = 1
        
        while restart_condition > 0:
            
            grompp_command = 'gmx grompp -f ' + args.input + ' -c ' + \
                gro_file  + ' -p ' + topology_file + ' -o ' + tpr_file
            mdrun_command = 'gmx mdrun -s ' + tpr_file + \
                ' -v -x ' + xtc_file + ' -c ' + gro_file_2
            subprocess.call(grompp_command, shell=True)
            subprocess.call(mdrun_command, shell=True)
            
            # Checking if step files are generated, indicating system failure due to lincs warnings
            file_list = os.listdir()
            step_files = [fn for fn in file_list if 'step' in fn]
            nr_step_files = len(step_files)
            restart_condition = nr_step_files
            if restart_condition > 0:
                subprocess.call('rm *step*', shell=True)
       
        # The old .gro file is saved under a new (unique) name
        fm.name_generator(gro_file, shock)
                
        # =====================================================================
        # Identifying water particles to remove or replace
        # =====================================================================
        
        u = mda.Universe(gro_file_2)
        all_atoms = u.select_atoms('all')
        water_atoms = u.select_atoms(water_name)
        water_pos = water_atoms.positions
        lipids = u.select_atoms('not resname W')
        nr_lipid_atoms = len(lipids)
        tj = u.trajectory[0]
        box_dimensions = tj.dimensions
        bin_size_a = args.bins*10
                
        # Cluster class used to convert the system in bins and find bin water clusters
        cluster_box = Cluster(lipids, box_dimensions, bin_size_a, args.lipids)
        
        target_lipids = cluster_box.l_particle_selector()        
        target_lipids_pos = target_lipids.positions
        lipid_bins = cluster_box.bin_converter_l(target_lipids_pos)
        lipid_bins_single = cl.multiples_remover(lipid_bins)
        bin_system = cluster_box.bin_system_maker(lipid_bins_single)
        bin_clusters = cluster_box.cluster_finder(bin_system)
        
        # using class Remover the correct cluster is identified and the indices
        # of the particles to be removed are obtained
        identify = Identifier(bin_clusters, box_dimensions, bin_size_a, args.removed)
        
        water_clusters = identify.cluster_selecter()
        inner_cluster = water_clusters[0]
        outer_cluster = water_clusters[1]
        np.random.shuffle(inner_cluster)
        all_water_bins = identify.bin_converter_w(water_pos)
        if args.shocktype == 'hypertonic':
            remove_water = identify.index_finder(all_water_bins, inner_cluster)
            remove_water_global = wi.local_to_global(nr_lipid_atoms, remove_water)
        else:
            remove_water = identify.index_finder(all_water_bins, outer_cluster)
            remove_water_global = wi.local_to_global(nr_lipid_atoms, remove_water)
        
        # =====================================================================
        # Modifying input files to create the new situation
        # =====================================================================
        if args.shocktype == 'hypertonic':
            
            if args.replace == 'no':
                file_collection = Remover(topology_file, all_atoms, args.removed)
            
                file_collection.water_remover_gro(remove_water_global, gro_file)
                file_collection.water_remover_top(topology_file_2)
                
            else:
                file_collection = Mover(outer_cluster, all_atoms, bin_size_a, args.removed)
                
                replacement_bins = file_collection.replacement_bin_identifier()
                new_positions = file_collection.position_generator(replacement_bins)
                file_collection.water_replacement_gro(remove_water_global, new_positions, gro_file)
        else:
            file_collection = Mover(inner_cluster, all_atoms)
            
            replacement_bins = file_collection.replacement_bin_identifier(args.removed)
            new_positions = wrep.position_generator(replacement_bins, bin_size_a)
            file_collection.water_replacement_gro(remove_water_global, new_positions, gro_file)
        
        # =====================================================================
        # Renaming, saving and deleting files
        # =====================================================================
                
        # renaming and saving files
        if args.shocktype == 'hypertonic':
            if args.replace == 'no':
                fm.name_generator(topology_file, shock)
                rename_command = 'mv ' + topology_file_2 + ' ' + topology_file
                subprocess.call(rename_command, shell=True)
        fm.name_generator(tpr_file, shock)
                               
        # making and saving an .xtc file
        fm.xtc_maker(lipids, shock)
        subprocess.call('mv *.xtc shockfiles', shell=True)
                
        # Save one gro file of the vesicle only
        if shock == 0:
            keep_atoms = u.select_atoms('not resname W')
            keep_atoms.write('vesicle_lipids.gro')
            subprocess.call('mv vesicle_lipids.gro shockfiles', shell=True)
        
        # Removing all backups and log files
        subprocess.call('rm *#', shell=True)
        subprocess.call('rm md.log', shell=True)
            
        # =====================================================================
        # Vesicle data acquisition
        # =====================================================================
                    
        # Calculating the volume and surface area and writing to file
        if args.vesdata == 'yes':
            if shock ==0:
                with open('vesicle_data.txt', 'a') as data:
                    data.write('{0:10} {1:20} {2:20} {3:20}'.format('Shock nr', 
                        'inner area (nm^2)', 'outer area (nm^2)', 'volume (nm^3)'))
                    data.write('\n')
                    data.write('-'*70)
                    data.write('\n')
                    
            selection = 'name ' + args.cloudpoints
            broken_system = VolArea(lipids, box_dimensions, selection)
            vol_area = broken_system.volume_area()
            
            with open('vesicle_data.txt', 'a') as data:
                
                data.write('{0:<10d} {1:<20.2f} {2:<20.2f} {3:<20.2f}'.format(shock, 
                    vol_area[0]/100, vol_area[1]/100, vol_area[2]/1000))
                data.write('\n')
                
    # =========================================================================
    # Post processing
    # =========================================================================
                  
    # Concatenating the xtc files
    concatcommand = 'gmx trjcat -cat -f shockfiles/*.xtc -o shockfiles/vesicle_lipids.xtc'
    subprocess.call(concatcommand, shell=True)
    subprocess.call('rm shockfiles/*t.xtc', shell=True)
  
if __name__ == '__main__':
    main()
