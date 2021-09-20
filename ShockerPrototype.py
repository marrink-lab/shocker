
import MDAnalysis as mda
import subprocess
import argparse

parser = argparse.ArgumentParser(description='Osmotic shock simulator')
parser.add_argument('-in','--input', help='mdp file with input data')
parser.add_argument('-f','--coordinates', help='gro file')
parser.add_argument('-t','--topology', help='topology file')
parser.add_argument('-r','--removed', type=int, help='nr of removed particles per shock')
parser.add_argument('-nr','--number', type=int, help='total number of shocks')
args = parser.parse_args() 

# %%

def water_selecter(gro_file):
    
    u = mda.Universe(gro_file)
    lipid_atoms = u.select_atoms('not resname W')
    water_atoms = u.select_atoms('resname W')
    nr_lipid_atoms = len(lipid_atoms)
    nr_water_atoms = len (water_atoms)
    
    # Calculating the center of mass (z-coordinate) of both bilayers
    top_atoms = lipid_atoms[int(nr_lipid_atoms/2):int(nr_lipid_atoms)]
    bottom_atoms = lipid_atoms[0:int(nr_lipid_atoms/2)]

    upper_limit_z = top_atoms.centroid()[2]
    lower_limit_z = bottom_atoms.centroid()[2]
    
    # Finding the indices of the water particles situated between the two bilayers
    water_inside = []

    for i in range(nr_water_atoms):
        water_z_pos = water_atoms.positions[i][2]
        if lower_limit_z < water_z_pos < upper_limit_z:
            water_inside.append(i)

      
    return water_inside

# %%

def water_remover_gro(old_gro, new_gro, positions):
    
    u = mda.Universe(old_gro)
    removed_atoms_li = positions
    water_atoms = u.select_atoms('resname W')
    all_atoms = u.select_atoms('all')
    nr_atoms = len(all_atoms)
    nr_water_atoms = len (water_atoms)
    no_water = nr_atoms - nr_water_atoms
    
    # Creating list of global indices of the atoms that have to be removed
    removed_atoms_gi = []                             
    for i in removed_atoms_li:
        removed_atoms_gi.append(i + no_water)
    removed_atoms_gi.append(-1)
    
    # Creating a new .gro file using the old .gro file and the positions of particles that have to be removed
    keep_atoms = u.select_atoms('index 0')
    for i in range(1,nr_atoms):
        if i != removed_atoms_gi[0]:
            number = 'index' + ' ' + str(i)
            keep_atoms = keep_atoms + u.select_atoms(number)
        else:
            removed_atoms_gi.pop(0)

    keep_atoms.write(new_gro)

# %%

def water_remover_top(old, new, nr_removed):

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
    no_ext = old.split('.')[0]
    new_top_name = no_ext + '_' + 's' + str(shock_nr) + '.top'

    return new_top_name

# %%

def gro_name_generator(old, shock_nr):
    no_ext = old.split('.')[0]
    new_gro_name = no_ext + '_' + 's' + str(shock_nr) + '.gro'
    
    return new_gro_name 

# %%

def extension_remover(filename):
    sp = filename.split('.')
    no_ext = sp[0]
    return no_ext

# %%

def main():
    nr_removed_atoms = args.removed
    
    for shock in range(args.number):
        
        input_file = args.input
        gro_file = args.coordinates
        tpr_file = extension_remover(gro_file) + '.tpr'
        gro_file_2 = extension_remover(gro_file) + '2.gro'
        xtc_file = extension_remover(gro_file) + '.xtc'
        topology_file = args.topology
        topology_file_2 = extension_remover(topology_file) + '2.top'
        
        # Executing the actual simulation     
        grompp_command = 'gmx grompp -f ' + input_file + ' -c ' + gro_file  + ' -p ' + topology_file + ' -o ' + tpr_file
        mdrun_command = 'gmx mdrun -s ' + tpr_file + ' -v -x ' + xtc_file + ' -c ' + gro_file_2
        subprocess.call(grompp_command, shell=True)
        subprocess.call(mdrun_command, shell=True)
        
        # Remove old .gro file
        remove_gro_command = 'rm ' + gro_file
        subprocess.call(remove_gro_command, shell=True)
        
        # Identifying the water particles inside
        inner_water = water_selecter(gro_file_2)
        
        # Removing water particles from the .gro file and the .top file
        removed_atoms = inner_water[0:nr_removed_atoms]
        water_remover_gro(gro_file_2, gro_file, removed_atoms)
        water_remover_top(topology_file, topology_file_2, nr_removed_atoms)
        
        # The old topology file is saved under a new (unique) name
        top_name = top_name_generator(topology_file, shock)
        t_name_change_command = 'mv ' + topology_file + ' ' + top_name
        subprocess.call(t_name_change_command, shell=True)
        
        # The new topology file is renamed to be used in the next cycle
        rename_command = 'mv ' + topology_file_2 + ' ' + topology_file
        subprocess.call(rename_command, shell=True)
        
        # We do not need the .xtc file
        remove_xtc_command = 'rm ' + xtc_file
        subprocess.call(remove_xtc_command, shell=True)
        
        # The old .gro file is saved under a new (unique) name
        gro_name = gro_name_generator(gro_file_2, shock)
        g_name_change_command = 'mv ' + gro_file_2 + ' ' + gro_name
        subprocess.call(g_name_change_command, shell=True)
        
        
if __name__ == '__main__':
    main()