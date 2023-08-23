# Shocker
Shocker is a python package for simulating the effects of water efflux (influx) during a hypertonic (hypotonic) osmotic shock by relocating water particles from the inner to the outer compartment (or vice versa). The algortihm can be applied to any structure containing an enclosed compartment of solvent such as a vesicle, tube or double bilayer system. The tool works with the md code GROMACS and can be applied to both AA systems and CG systems. Shocker is applicable to membranes containing proteins and other molecules such as polymers and DNA. Shocker divides the simulation in pumping cycles, defined as relocation of a number of water particles followed by relaxation time.

# Instructions

## Installation
Clone repository
```
git clone https://github.com/marrink-lab/shocker
```

Create virtual environment

```python3 -m venv venv```

On Peregrine and Snellius load python module

Peregrine

```module load Python/3.8.6-GCCcore-10.2.0```

Snellius
```
module load 2021
module load Python/3.9.5-GCCcore-10.3.0
```

Install from setup folder

```pip install ./```

In bash script add:

```
module load 2021
module load Python/3.9.5-GCCcore-10.3.0
source ~/venv/bin/activate
```

## Required files
In order to perform an osmotic shock simulation Shocker needs an equilibrated simulation box (.GRO) containing a structure with a closed water compartment (enclosed by membrane or periodic boundaries). Besides, a parameter file (.MDP), a topology file(.TOP) and an index file (.NDX) should be provided, where the index file should discern membrane, solvent and ion (if present) particle groups.

## Example 1: Hypertonic shock on a vesicle
A simulation of a hypertonic osmotic shock on a vesicle is executed with the following command:
```
shocker -c POPC_vesicle.gro -f martini_md_anis.mdp -p topol_POPC.top -r 200 -e 200 -n index_POPC.ndx
```
(-r) indicates the number of solvent particles relocated each pumping cycle  
(-e) is the last cycle

This is the minimal information Shocker needs to perform a simulation. The default relaxation time after each pumping cycle is 2 ns (100.000 timesteps with a size of 0.02 ps). The stepsize and the number of steps can be modified with the (-pts) flag and the (-it) flag respectively. The default names of the structure groups are 'membrane', 'water' and 'ions'. These names
should be used in the index file as well. It is possible to use other names for these groups, but this should be provided with the (-mem), (-water) and (-ions) flags.
When the process is finished the simulation directory contains a folder called 'shockfiles' where the .tpr and .gro files of each pumping cycle are stored, as well as the trajectory of the complete simulation.

## Example 2: Shape analysis during hypertonic shock on a vesicle
The Shocker package contains a few analysis methods that enables the user to monitor the shape of a vesicle over the course of the osmotic shock simulation. In this case we consider a simulation box containing a vesicle in a solvent with ions. The simulation and data calculation is executed with the following command:
```
shocker -c POPC_ion.gro -f martini_md_anis.mdp -p topol_ion.top -r 200 -e 200 -n index_ion.ndx -vd yes -i yes
```
(-vd) indicates if vesicle data should be calculated  
(-i) indicates the presence of ions  

If both of these values are set to 'yes' a file named 'vesicle_data.txt' will be generated containing in the first place the vesicle's inner area, outer area and volume. Besides, the reduced area difference, reduced volume and a parameter indicating the sphericity will be calculated. In case ions are present the concentration inside and outside the vesicle will be calculated.

## Example 3: Hypotonic osmotic shock
Except for deflating vesicles during a hypertonic shock, Shocker is capable of performing hypotonic shock simulations to study the effects of volume increase. Such a simulation can be initiated as follows:
```
shocker -c POPC_ion.gro -f martini_md_anis.mdp -p topol_ion.top -r 200 -e 200 -n index_ion.ndx -shock hypotonic
```

## Restart/continue
When initiating a new simulation Shocker always starts with iteration 0. If a pumping simulation is executed in a folder that already contains pumping data, Shocker continues from the last pumping cycle (e.g. in case of sudden termination).

## Extra mdrun options
As default the mdrun command is executed using the minimum amount of imput options. Especially when running Shocker on a cluster the user may want to add some extra options. This can be done by providing an option file (text file) with the -options flag, for example:

```
-ntmpi 6
-ntomp 3
-dlb auto
-pin on
-nb gpu
-bonded gpu
```

