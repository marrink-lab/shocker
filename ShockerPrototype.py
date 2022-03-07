
import MDAnalysis as mda
import subprocess
import argparse
import numpy as np
import MDAnalysis.analysis.leaflet
import os
import pyvista as pv

parser = argparse.ArgumentParser(description='Osmotic shock simulator', add_help=True)
parser.add_argument('-in','--input', help='mdp file with input data')
parser.add_argument('-f','--coordinates', help='gro file')
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
                    
def mdp_value_changer(old, nr_its, valnm):
    '''Changes value in mdp file'''
    with open(old) as cur_mdp:
            lines_mdp = cur_mdp.readlines()
            nr_mdp_lines = len(lines_mdp)
            for i in range(nr_mdp_lines):
                spl = lines_mdp[i].split(' ')
                if spl[0] != valnm:
                    with open('new.mdp', 'a') as new_mdp:
                        new_mdp.write(lines_mdp[i])
                else:
                    spl[-1] = str(nr_its) + '\n'
                    new_spl = ' '.join(spl)
                    with open('new.mdp', 'a') as new_mdp:
                        new_mdp.write(new_spl)

    mdp_remove_command = 'rm ' + old
    mdp_name_change_command = 'mv ' + 'new.mdp ' + old                   

    subprocess.call(mdp_remove_command, shell=True)                   
    subprocess.call(mdp_name_change_command, shell=True)

# %%

def name_generator(old, shock_nr, ext):
    '''Generates a name to save the top file every cycle'''
    no_ext = old.split('.')[0]
    new_name = no_ext + '_' + 's' + str(shock_nr) + ext

    return new_name

# %%

def extension_remover(filename):
    ''''Removes the extension of a filename'''
    sp = filename.split('.')
    no_ext = sp[0]
    return no_ext

# %%
    
def l_particle_selector(lipidlist):
    '''selects the lipid particles targeted for creating the lipid bins according to the lipid names'''
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

def center_of_geometry(points):
    '''calculates the center of geometry'''
    nr = len(points)
    x = 0
    y = 0
    z = 0
    
    for i in points:
        x = x + i[0]
        y = y + i[1]
        z = z + i[2]
    
    com = np.array([x/nr, y/nr, z/nr])
    
    return com

# %%

def boundary_corrector(lips, box_dimensions):
    '''Correction method for structures passing the periodic boundary'''
    l = MDAnalysis.analysis.leaflet.LeafletFinder(lips, 'name GL1')
    
    group_len = []
    for i in range(len(l.groups())):
        group_len.append(len(l.group(i)))
    list.sort(group_len)
    
    base_index = []
    for i in range(len(l.groups())):
        if len(l.group(i)) >= group_len[-2]:
            base_index.append(i)
    base = l.group(base_index[0]) + l.group(base_index[1])
    
    for i in l.groups():
        if len(i) < group_len[-2]:
            
            cluster = i.positions
            
            com = center_of_geometry(cluster)
            
            absvalues = list(abs(box_dimensions[:3] - com)) + list(abs([0,0,0] - com))
            
            minimum = np.min(absvalues)
            
            min_index = list(absvalues).index(minimum)
    
            if min_index == 0:
                for v in cluster:
                    v[0] = v[0] - box_dimensions[0]
    
            if min_index == 1:
                for v in cluster:
                    v[1] = v[1] - box_dimensions[1]
    
            if min_index == 2:
                for v in cluster:
                    v[2] = v[2] - box_dimensions[2]
    
            if min_index == 3:
                for v in cluster:
                    v[0] = v[0] + box_dimensions[0]
    
            if min_index == 4:
                for v in cluster:
                    v[1] = v[1] + box_dimensions[1]
    
            if min_index == 5:
                for v in cluster:
                    v[2] = v[2] + box_dimensions[2]  

            i.positions = cluster
            
            base = base + i
    return base

# %%

def mesh_maker(cloud):
    '''Generates a triangulated mesh based on a pointcloud'''
    points = pv.wrap(cloud)
    surf = points.reconstruct_surface(nbr_sz=6)
    
    return surf

# %%

def volume_area(leafl, pointc, boxdim):
    '''Calculates the area of both leaflets and the volume of the vesicle'''
    leaflets = leafl
    pointcloud = pointc
    outer = leaflets.groups(0).positions
    inner = leaflets.groups(1).positions
    
    if len(leaflets.groups()) > 2:
                    
        pointcloud = boundary_corrector(pointcloud, boxdim)
        leaflets = MDAnalysis.analysis.leaflet.LeafletFinder(pointcloud, 'name GL1')
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

# %%

def xtc_maker(lipids, cycle):
    '''Generates an xtc file from a gro file in order to concatenate the cycles'''
    new_gro_name = 'only_lipids_temp.gro'
    
    lipids.write(new_gro_name)
    
    if len(str(cycle)) == 1:
        new_name = 'vesicle_sA' + str(cycle) + '_t.xtc'
    elif len(str(cycle)) == 2:
        new_name = 'vesicle_sB' + str(cycle) + '_t.xtc'
    else:
        new_name = 'vesicle_sC' + str(cycle) + '_t.xtc'
    
    convcommand = 'gmx trjconv -f ' + new_gro_name + ' -o ' + new_name  
    subprocess.call(convcommand, shell=True)
    
    rmcommand = 'rm ' + new_gro_name
    subprocess.call(rmcommand, shell=True)
    
    return new_name

# %%

def main():
    nr_removed_atoms = args.removed
    
    # changing the number of iterations per cycle and frequency of position data storage (default is 100.000 iterations)
    mdp_value_changer(args.input, args.iterations, 'nsteps')
    mdp_value_changer(args.input, args.xtc, 'nstxout-compressed')
    
    input_file = args.input
    subprocess.call('mkdir shockfiles', shell=True)
    
    #changing the nstxout value (default is 100.000)
    
    for shock in range(args.start, args.end):
        
        gro_file = args.coordinates
        tpr_file = extension_remover(gro_file) + '.tpr'
        gro_file_2 = extension_remover(gro_file) + '2.gro'
        xtc_file = extension_remover(gro_file) + '.xtc'
        topology_file = args.topology
        topology_file_2 = extension_remover(topology_file) + '2.top'
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
        gro_name = name_generator(gro_file, shock, '.gro')
        g_name_change_command = 'mv ' + gro_file + ' ' + gro_name
        subprocess.call(g_name_change_command, shell=True)
        dir_change = 'mv ' + gro_name + ' shockfiles'
        subprocess.call(dir_change, shell=True)
        
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
        only_lipids = u.select_atoms('not resname W')
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
        top_name = name_generator(topology_file, shock, '.top')
        t_name_change_command = 'mv ' + topology_file + ' ' + top_name
        subprocess.call(t_name_change_command, shell=True)
        dir_change = 'mv ' + top_name + ' shockfiles'
        subprocess.call(dir_change, shell=True)
        
        # The new topology file is renamed to be used in the next cycle
        rename_command = 'mv ' + topology_file_2 + ' ' + topology_file
        subprocess.call(rename_command, shell=True)
        
        # Removing log file
        subprocess.call('rm md.log', shell=True)
        
        # making and saving an .xtc file
        xtc_file = xtc_maker(only_lipids, shock)
        dir_change = 'mv ' + xtc_file + ' shockfiles'
        subprocess.call(dir_change, shell=True)
        
        # Save one gro file of the vesicle only
        if shock == 0:
            keep_atoms = u.select_atoms('not resname W')
            keep_atoms.write('vesicle_lipids.gro')
            subprocess.call('mv vesicle_lipids.gro shockfiles', shell=True)
                    
        # Saving the .tpr file
        tpr_name = name_generator(tpr_file, shock, '.tpr')
        tpr_name_change_command = 'mv ' + tpr_file + ' ' + tpr_name
        subprocess.call(tpr_name_change_command, shell=True)
        dir_change = 'mv ' + tpr_name + ' shockfiles'
        subprocess.call(dir_change, shell=True)
        
        # Removing all other backups
        subprocess.call('rm *#', shell=True)
        
        # Calculating the volume and surface area and writing to file
        if args.vesdata == 'yes':
            selection = 'name ' + args.cloudpoints
            leaflets = MDAnalysis.analysis.leaflet.LeafletFinder(u, selection)
            pointcloud = u.select_atoms(selection).positions
            vol_area = volume_area(leaflets, pointcloud, box_dimensions)
            
            with open('vesicle_data.txt', 'a') as data:
                
                data.write(str(shock))
                data.write(' ')
                data.write(str(vol_area[0]))
                data.write(' ')
                data.write(str(vol_area[1]))
                data.write(' ')
                data.write(str(vol_area[2]))
                data.write('\n')
    
    # Concatenating the xtc files
    concatcommand = 'gmx trjcat -cat -f shockfiles/*.xtc -o shockfiles/vesicle_lipids.xtc'
    subprocess.call(concatcommand, shell=True)
    subprocess.call('rm shockfiles/*t.xtc', shell=True)
    
if __name__ == '__main__':
    main()
    
    
    
    
    