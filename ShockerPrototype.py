
import MDAnalysis as mda
import subprocess
import argparse
import numpy as np
import MDAnalysis.analysis.leaflet

parser = argparse.ArgumentParser(description='Osmotic shock simulator')
parser.add_argument('-in','--input', help='mdp file with input data')
parser.add_argument('-f','--coordinates', help='gro file')
parser.add_argument('-t','--topology', help='topology file')
parser.add_argument('-r','--removed', type=int, help='nr of removed particles per shock')
parser.add_argument('-e','--end', type=int, help='nr of the last cycle')
parser.add_argument('-w','--water', default='resname W', help='name of water particles')
parser.add_argument('-b','--bins', type=int, default=0, help='increase or decrease the bin size (standard size is 1nm^3)')
parser.add_argument('-s','--start', type=int, default=0, help='define the starting cycle')
parser.add_argument('-x','--xtc', default='delete', help='save (save) xtc or not (delete)')
parser.add_argument('-lip','--lipids', nargs='+', help='the names of the lipids the membrane consists of')
args = parser.parse_args() 

# %%

def bin_converter(pos, boxdim, nrbins):
    '''Assigns each coordinate to a bin of a chosen size'''
    x = np.ceil(pos[:,0]/(boxdim[0]/nrbins[0]))
    y = np.ceil(pos[:,1]/(boxdim[1]/nrbins[1]))
    z = np.ceil(pos[:,2]/(boxdim[2]/nrbins[2]))
    
    binsfloat = np.stack((x,y,z), axis=-1)
    binsint = np.int_(binsfloat)
    
    for i in binsint:
        for u in range(3):
            if i[u] > nrbins[u]:
                i[u] = i[u] - nrbins[u]
            if i[u] < 0:
                i[u] = nrbins[u] + i[u]
       
    return binsint

# %%
    
def multiples_remover(allbins):
    '''Makes sure that each value appears only once in the list'''
    bintuples = [tuple(row) for row in allbins]
    binuniques = np.unique(bintuples, axis=0)
    
    return binuniques

# %%
    
def bin_system_maker(lbins, nrbins):
    '''A bin system is created in which 0 is a water bin and 1 is a lipid bin'''
    bsystem = np.zeros((nrbins[0], nrbins[1], nrbins[2]))
    
    for i in lbins:
        bsystem[i[0]-1][i[1]-1][i[2]-1] = 1
    
     
    return bsystem

# %%

def zero_finder(system):
    ''''Finds the position of a zero in a cubic bin system. 0 is returned if no zero is found'''
    zero = 0
    for x in range(len(system)):
        for y in range(len(system[x])):
            for z in range(len(system[x][y])):
                if system[x][y][z] == 0:
                    zero = (x,y,z)
    return zero

# %%

def index_finder(w_all, w_clstr, nrremoved):
    '''Of each element in w_clstr it returns the index in w_all, until enough particles have been collected'''
    indices  = []
    i = 0
    while len(indices) < nrremoved:
        bin_index = np.where((w_all[:,0] == w_clstr[i][0]) & (w_all[:,1] == w_clstr[i][1]) & (w_all[:,2] == w_clstr[i][2]))[0] 
        for u in bin_index:
            indices.append(int(u))
        i = i + 1
        
    return indices[:nrremoved]

# %%
    
def cluster_finder(binsystem, bindim):
    '''Çlusters the zeros in binsystem and returns all found clusters'''
    
    binx = bindim[0]
    biny = bindim[1]
    binz = bindim[2]
    
    black_list = binsystem
    cluster_list = []
    
    while zero_finder(black_list) != 0:
        
        next_start = zero_finder(black_list)
        
        gray_list = [next_start]
        white_list = []
        temp_list = []
        
        while len(gray_list) > 0:
            for u in gray_list:
                if u[0] < (binx - 1):
                    
                    if black_list[u[0] + 1][u[1]][u[2]] == 0:
                        temp_list.append((u[0] + 1, u[1], u[2]))
                        black_list[u[0] + 1][u[1]][u[2]] = -1
                        
                if u[0] > 0:
                    
                    if black_list[u[0] - 1][u[1]][u[2]] == 0:
                        temp_list.append((u[0] - 1, u[1], u[2]))
                        black_list[u[0] - 1][u[1]][u[2]] = -1
                        
                if u[1] < (biny - 1):
                    
                    if black_list[u[0]][u[1] + 1][u[2]] == 0:
                        temp_list.append((u[0], u[1] + 1, u[2]))
                        black_list[u[0]][u[1] + 1][u[2]] = -1
                        
                if u[1] > 0:
                    
                    if black_list[u[0]][u[1] - 1][u[2]] == 0:
                        temp_list.append((u[0], u[1] - 1, u[2]))
                        black_list[u[0]][u[1] - 1][u[2]] = -1
                        
                if u[2] < (binz - 1):
                    
                    if black_list[u[0]][u[1]][u[2] + 1] == 0:
                        temp_list.append((u[0], u[1], u[2] + 1))
                        black_list[u[0]][u[1]][u[2] + 1] = -1
                        
                if u[2] > 0:
                    
                    if black_list[u[0]][u[1]][u[2] - 1] == 0:
                        temp_list.append((u[0], u[1], u[2] - 1))
                        black_list[u[0]][u[1]][u[2] - 1] = -1
                white_list.append(u)
                black_list[u[0]][u[1]][u[2]] = -1
                
            gray_list = temp_list
            temp_list = []
        cluster_list.append(white_list)
    
    
    return cluster_list

# %%

def cluster_selecter(clist, nrbins):
    ''''Selects the relevant cluster'''
    
    target_cluster = []
    for i in range(len(clist)):
        if (0,0,0) not in clist[i] and (nrbins[0]-1, nrbins[1]-1, nrbins[2]-1) not in clist[i]:
            target_cluster = target_cluster + clist[i]
    
            
    return target_cluster

# %%
    
def water_selecter(lipids, water, nrbins, boxdim, nrremoved):
    '''Selects the incides of the water particles that have to be removed'''
    lipid_bins = bin_converter(lipids, boxdim, nrbins)
    
    lipid_bins_single = multiples_remover(lipid_bins)
    
    bin_system = bin_system_maker(lipid_bins_single, nrbins)
    
    clusters = cluster_finder(bin_system, nrbins)
    
    inner_water_cluster = cluster_selecter(clusters, nrbins)
    
    np.random.shuffle(inner_water_cluster)
    
    water_bins = bin_converter(water, boxdim, nrbins)
    
    removed = index_finder(water_bins, inner_water_cluster, nrremoved)
    
    return removed

# %%

def local_to_global(no_water, localindex):
    ''''(local) indices found in the list of water particles are converted to global indices of the complete system list'''
    global_index = []
    for i in localindex:
        global_index.append(i + no_water)
        
    return global_index

# %%

def not_string(positions):
    '''A string is created of particle indices that have to be removed from the system'''
    nstring = 'not index '
    for i in positions:
        nstring = nstring + str(i) + ' '
        
    return nstring

# %%

def water_remover_top(old, new, nr_removed):
    '''Çreates a new topology file with the new number of water particles'''
    # The original topology file is modified to make sure it shows the correct number of W-particles
    with open(old) as cur_top:
        lines_top = cur_top.readlines()
        nr_top_lines = len(lines_top)
            
        # All lines are copied to the new .top file except for the one containing the number of W-particles.
        for i in range(nr_top_lines):
            sp_string_top = lines_top[i].split(' ')
            if sp_string_top[0] != 'W':
                with open(new, 'a') as new_top:
                    new_top.write(lines_top[i])
            else:
                  
                # The number of removed W-particles has to be subtracted from the current number.
                cur_nr_water = sp_string_top[1].split('\n')
                new_nr_water = int(cur_nr_water[0]) - nr_removed
                new_string_element_water = str(new_nr_water) + '\n'
            
                sp_string_top[1] = new_string_element_water
                new_string_water = ' '.join(sp_string_top)
            
                with open(new, 'a') as new_top:
                    new_top.write(new_string_water)

# %%

def top_name_generator(old, shock_nr):
    '''Generates a name to save the top file every cycle'''
    no_ext = old.split('.')[0]
    new_top_name = no_ext + '_' + 's' + str(shock_nr) + '.top'

    return new_top_name

# %%

def gro_name_generator(old, shock_nr):
    '''Generates a name to save the gro file every cycle'''
    no_ext = old.split('.')[0]
    new_gro_name = no_ext + '_' + 's' + str(shock_nr) + '.gro'
    
    return new_gro_name

# %%

def xtc_name_generator(old, shock_nr):
    '''Generates a name to save the xtc file every cycle'''
    no_ext = old.split('.')[0]
    new_xtc_name = no_ext + '_' + 's' + str(shock_nr) + '.xtc'
    
    return new_xtc_name 

# %%

def tpr_name_generator(old, shock_nr):
    '''Generates a name to save the tpr file every cycle'''
    no_ext = old.split('.')[0]
    new_tpr_name = no_ext + '_' + 's' + str(shock_nr) + '.tpr'
    
    return new_tpr_name 

# %%

def extension_remover(filename):
    ''''Removes the extension of a filename'''
    sp = filename.split('.')
    no_ext = sp[0]
    return no_ext

# %%
    
def l_particle_selector(lipidlist):
    
    lipids = lipidlist
    tail_particles = []
    
    for lipid in lipids:
        
        with open('Martini3.LIB', 'r') as lipfile:
            lines = lipfile.readlines()
            for i in range(len(lines)):
                name = lines[i].strip()[1:-1]
                
                if lipid == 'CHOL':
                    tail_particles.append('C1 C2')
                
                elif name.replace(' ','') == lipid:
                                        
                    lipid_particles = []
                    c = 1
                    while lines[i+c] != '\n':
                        lipid_particles.append(lines[i+c].split()[1])
                        c = c + 1
                    
                    tail_length = int(lipid_particles[-1][1])
                    tail_positions = [tail_length - 2, tail_length - 1, tail_length]
                    
                    for atom in lipid_particles:
                        if atom[1] == str(tail_positions[0]) or atom[1] == str(tail_positions[1]) or atom[1] == str(tail_positions[2]):
                            tail_particles.append(atom)
                    
    unique_tail_particles = []
    for i in tail_particles:
        if i not in unique_tail_particles:
            unique_tail_particles.append(i)
    
    particle_string = ''
    for i in unique_tail_particles:
        particle_string = particle_string + ' ' + str(i)
        
    selection_command = 'name' + particle_string  
    
    return selection_command

# %%

def main():
    nr_removed_atoms = args.removed
    
    for shock in range(args.start, args.end):
        
        input_file = args.input
        gro_file = args.coordinates
        tpr_file = extension_remover(gro_file) + '.tpr'
        gro_file_2 = extension_remover(gro_file) + '2.gro'
        xtc_file = extension_remover(gro_file) + '.xtc'
        topology_file = args.topology
        topology_file_2 = extension_remover(topology_file) + '2.top'
        water_name = args.water
        
        # Executing the actual simulation     
        grompp_command = 'gmx grompp -f ' + input_file + ' -c ' + gro_file  + ' -p ' + topology_file + ' -o ' + tpr_file
        mdrun_command = 'gmx mdrun -s ' + tpr_file + ' -v -x ' + xtc_file + ' -c ' + gro_file_2
        subprocess.call(grompp_command, shell=True)
        subprocess.call(mdrun_command, shell=True)
        
        # The old .gro file is saved under a new (unique) name
        gro_name = gro_name_generator(gro_file, shock)
        g_name_change_command = 'mv ' + gro_file + ' ' + gro_name
        subprocess.call(g_name_change_command, shell=True)
        
        # Variables needed for the calculations
        u = mda.Universe(gro_file_2)
        water_atoms = u.select_atoms(water_name)
        water_pos = water_atoms.positions
        
        # Lipid particles used to define the lipid bins are selected
        membrane_composition = args.lipids
        lip_part_selection = l_particle_selector(membrane_composition)
        target_lipids = u.select_atoms(lip_part_selection)
        
        target_lipids_pos = target_lipids.positions
        tj = u.trajectory[0]
        box_dimensions = tj.dimensions
        all_atoms = u.select_atoms('all')
        nr_atoms = len(all_atoms)
        nr_water_atoms = len(water_atoms)
        no_water = nr_atoms - nr_water_atoms
        bin_multiplier = (100-args.bins)/1000
        nr_bins = [int(box_dimensions[0]*bin_multiplier),int(box_dimensions[1]*bin_multiplier),int(box_dimensions[2]*bin_multiplier)]
        
        # Identifying the indices of the water particles inside we want to remove
        removed_atoms = water_selecter(target_lipids_pos, water_pos, nr_bins, box_dimensions, nr_removed_atoms)
        global_indices = local_to_global(no_water, removed_atoms)
        
        # Removing water particles from the .gro file and the .top file
        remove_string = not_string(global_indices)
        keep_atoms = u.select_atoms(remove_string)
        keep_atoms.write(gro_file)
        
        water_remover_top(topology_file, topology_file_2, nr_removed_atoms)
        
        # The old topology file is saved under a new (unique) name
        top_name = top_name_generator(topology_file, shock)
        t_name_change_command = 'mv ' + topology_file + ' ' + top_name
        subprocess.call(t_name_change_command, shell=True)
        
        # The new topology file is renamed to be used in the next cycle
        rename_command = 'mv ' + topology_file_2 + ' ' + topology_file
        subprocess.call(rename_command, shell=True)
        
        # Removing log file
        subprocess.call('rm md.log', shell=True)
        
        # Saving or removing .xtc file
        if args.xtc == 'save':
            xtc_name = xtc_name_generator(xtc_file, shock)
            x_name_change_command = 'mv ' + xtc_file + ' ' + xtc_name
            subprocess.call(x_name_change_command, shell=True)
        else:
            remove_xtc_command = 'rm ' + xtc_file
            subprocess.call(remove_xtc_command, shell=True)
            
        # Saving the .tpr file
        tpr_name = tpr_name_generator(tpr_file, shock)
        tpr_name_change_command = 'mv ' + tpr_file + ' ' + tpr_name
        subprocess.call(tpr_name_change_command, shell=True)
        
        # Removing all other backups
        subprocess.call('rm *#', shell=True)
        
if __name__ == '__main__':
    main()