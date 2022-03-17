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
import os
import Shocker.shocker.src.cluster as cl
import Shocker.shocker.src.file_mod as fm
import Shocker.shocker.src.vesicle_data as vd
import Shocker.shocker.src.water_remover as wr

def main():
    
    parser = argparse.ArgumentParser(description='Osmotic shock simulator', add_help=True)
    parser.add_argument('-in','--input', help='mdp file with input data')
    parser.add_argument('-f','--name', help='name of the files')
    parser.add_argument('-t','--topology', help='topology file')
    parser.add_argument('-r','--removed', type=int, help='nr of removed particles per shock')
    parser.add_argument('-e','--end', type=int, help='nr of the last cycle')
    parser.add_argument('-w','--water', default='resname W', help='name of water particles')
    parser.add_argument('-b','--bins', type=int, default=0, help='increase or decrease the bin size (standard size is 1nm^3)')
    parser.add_argument('-s','--start', type=int, default=0, help='define the starting cycle')
    parser.add_argument('-x','--xtc', type=int, default=0, help='frequency of storing position data (default = 0)')
    parser.add_argument('-lip','--lipids', nargs='+', help='the names of the lipids the membrane consists of')
    parser.add_argument('-it', '--iterations', type=int, default=100000, help='the number of iterations per cycle')
    parser.add_argument('-min','--enermin', default='no', help='energy minimization after each pumping cycle (yes or no)')
    parser.add_argument('-emnr', '--enernr', type=int, default=1, help='number of energy minimization rounds per cycle')
    parser.add_argument('-vd', '--vesdata', default='no', help='whether or not leaflet area and vesicle volume is calculated and stored each cycle')
    parser.add_argument('-cp', '--cloudpoints', default='GL1', help='defines which atoms are used to construct a triangulated surface of the vesicle')
    args = parser.parse_args()
    
    nr_removed_atoms = args.removed
    
    # changing the number of iterations per cycle and frequency of position data storage (default is 100.000 iterations)
    fm.mdp_value_changer(args.input, 'nsteps', args.iterations)
    fm.mdp_value_changer(args.input, 'nstxout-compressed', args.xtc)
    
    input_file = args.input
    subprocess.call('mkdir shockfiles', shell=True)
    
    #changing the nstxout value (default is 100.000)
    
    for shock in range(args.start, args.end):
        
        gro_file = args.name + '.gro'
        tpr_file = args.name + '.tpr'
        gro_file_2 = args.name + '2.gro'
        xtc_file = args.name + '.xtc'
        topology_file = args.topology + '.top'
        topology_file_2 = args.topology + '2.top'
        water_name = args.water
        
        # Energy minimization
        if args.enermin == 'yes':
            
            for i in range(args.enernr):
                grompp_command_em = 'gmx grompp -f min.mdp -c ' + gro_file  + ' -p ' + topology_file + ' -maxwarn 1'
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
            
            grompp_command = 'gmx grompp -f ' + input_file + ' -c ' + gro_file  + ' -p ' + topology_file + ' -o ' + tpr_file
            mdrun_command = 'gmx mdrun -s ' + tpr_file + ' -v -x ' + xtc_file + ' -c ' + gro_file_2
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
        wr.name_generator(gro_file, shock)
                
        # Variables needed for the calculations
        u = mda.Universe(gro_file_2)
        water_atoms = u.select_atoms(water_name)
        water_pos = water_atoms.positions
        
        # Lipid particles used to define the lipid bins are selected
        membrane_composition = args.lipids
        lip_part_selection = cl.l_particle_selector(membrane_composition)
        target_lipids = u.select_atoms(lip_part_selection)
        
        target_lipids_pos = target_lipids.positions
        tj = u.trajectory[0]
        box_dimensions = tj.dimensions
        all_atoms = u.select_atoms('all')
        nr_atoms = len(all_atoms)
        nr_water_atoms = len(water_atoms)
        no_water = nr_atoms - nr_water_atoms
        only_lipids = u.select_atoms('not resname W')
        bin_multiplier = (100-args.bins)/1000 #factor with which the box dimensions are multiplied to obtain the number of bins
        nr_bins = [int(box_dimensions[0]*bin_multiplier),int(box_dimensions[1]*bin_multiplier),int(box_dimensions[2]*bin_multiplier)]
        
        # Identifying the indices of the water particles inside we want to remove
        removed_atoms = wr.water_selecter(target_lipids_pos, water_pos, nr_bins, box_dimensions, nr_removed_atoms)
        global_indices = wr.local_to_global(no_water, removed_atoms)
        
        # Removing water particles from the .gro file and the .top file
        remove_string = wr.not_string(global_indices)
        keep_atoms = u.select_atoms(remove_string)
        keep_atoms.write(gro_file)
        
        wr.water_remover_top(topology_file, topology_file_2, nr_removed_atoms)
        
        # The old topology file is saved under a new (unique) name
        wr.name_generator(topology_file, shock)
                
        # The new topology file is renamed to be used in the next cycle
        rename_command = 'mv ' + topology_file_2 + ' ' + topology_file
        subprocess.call(rename_command, shell=True)
        
        # Removing log file
        subprocess.call('rm md.log', shell=True)
        
        # making and saving an .xtc file
        vd.xtc_maker(only_lipids, shock)
        subprocess.call('mv *.xtc shockfiles', shell=True)
                
        # Save one gro file of the vesicle only
        if shock == 0:
            keep_atoms = u.select_atoms('not resname W')
            keep_atoms.write('vesicle_lipids.gro')
            subprocess.call('mv vesicle_lipids.gro shockfiles', shell=True)
                    
        # Saving the .tpr file
        wr.name_generator(tpr_file, shock)
                
        # Removing all other backups
        subprocess.call('rm *#', shell=True)
        
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
            vol_area = vd.volume_area(only_lipids, box_dimensions, selection)
            
            with open('vesicle_data.txt', 'a') as data:
                
                data.write('{0:<10d} {1:<20.2f} {2:<20.2f} {3:<20.2f}'.format(shock, 
                    vol_area[0]/100, vol_area[1]/100, vol_area[2]/1000))
                data.write('\n')
                  
    # Concatenating the xtc files
    concatcommand = 'gmx trjcat -cat -f shockfiles/*.xtc -o shockfiles/vesicle_lipids.xtc'
    subprocess.call(concatcommand, shell=True)
    subprocess.call('rm shockfiles/*t.xtc', shell=True)
    
if __name__ == '__main__':
    main()
