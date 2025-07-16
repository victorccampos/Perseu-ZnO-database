import numpy as np
import os
from typing import Tuple

# Atomic Simulation Environment (ASE) imports
from ase.io import read                         # read ZnO_template.in
from ase import Atoms
from ase.build import make_supercell
from ase.io.espresso import write_espresso_in



def make_pw_input_from_existing_in(input_file_path: str, a: float, covera: float, cellsize: Tuple[int, int, int]):
    """
   Create a supercell input file for Quantum ESPRESSO based on an existing input file.
   Args:
        input_file_path (str): Path to the existing input file (e.g., ZnO_template.in).
        a (float): Lattice parameter 'a' in Angstroms.
        covera (float): c/a ratio for the wurtzite structure.
        cellsize (Tuple[int, int, int]): Tuple representing the supercell size in terms of unit cells (nx, ny, nz).
   Returns:
        None: The function creates a directory and writes the input file for the supercell. 
    """
    # OLD (wrong): bulk(name='ZnO', crystalstructure='wurtzite', a=a, covera=covera) 
    # --> creates convencional cell (hexagonal with 8 atoms)#

    primitive_cell = read(input_file_path, format='espresso-in') # Primitive wurtzite cell for ZnO from scripts/ZnO_template.in
    nx, ny, nz = cellsize
    supercell_matrix = np.array([
        [nx, 0, 0],
        [0, ny, 0],
        [0, 0, nz] ]
    )
    supercell = make_supercell(prim=primitive_cell, P=supercell_matrix)
    # Create the directory for the supercell type
    cell_dirname = f'cell{nx}{ny}{nz}.in'
    os.makedirs(cell_dirname, exist_ok=True)

    # Input file name based on the parameters
    input_name: str = f'ZnO-{a:.2f}-{covera:.2f}-{nx}{ny}{nz}.in'
    print(f'Produzindo : {input_name}')
    
    # Input comum para todas as simulações
    HEADER_INPUT = {
    'control': {
        'calculation': 'scf',
        'prefix': f'{input_name.replace(".in", "")}',
        'pseudo_dir': '/home/jvc/QEspresso7.2/ZnO_database/pseudos',
        'outdir': './',
        'disk_io': 'none',
        'verbosity': 'high',
        'tprnfor': True
    },
    'system': {
        'ibrav': 0,
        'ecutwfc': 50,
        'ecutrho': 400,
        'occupations': 'fixed'
    },
    'electrons': {
        'conv_thr': 1.0e-8,
        'mixing_beta': 0.3
    }
    }
    PSEUDOS = {'Zn': 'Zn.upf','O': 'O.upf'}
    
    # Adapt K-grid to larger supercells, mantaining grid density
    kpoints_original = np.array([6, 6, 6])  
    kpoints_supercell = kpoints_original // np.array(cellsize)  
    kpoints_supercell[kpoints_supercell == 0] = 1 # Ensure no zero k-points;

    # Write the input file using ASE's espresso interface
    write_espresso_in(
        #file = os.path.join(cell_dirname, input_name), # One folder per supercell
        file = input_name, # If you want to write in the current directory
        atoms = supercell,
        input_data = HEADER_INPUT,
        pseudopotentials = PSEUDOS,
        kpts = tuple(kpoints_supercell),
        koffset = (0,0,0),
        crystal_coordinates = True # TODO: check if this is needed
    )

# ============= Structural Parameters ============== #
percent_range = np.arange(start= -0.10, stop=0.12, step=0.02)

# vc-relax data from ZnO_system_preparation/data/vc-relax2
CELLDM1 = 6.178821408099141            # in Bohr units
CELLDM1_Angstroms = CELLDM1 * 0.52918  # in Angstroms
CELLDM3 = 1.614358356153010 # Adimensional c/a ratio

# Strained parameters - FUTURE
a_values = CELLDM1_Angstroms * (1 + percent_range)
c_values = CELLDM3 * (1 + percent_range)
# ================================================ #


if  __name__ == "__main__":
    # Path to the template input file - get 1x1x1 atoms structure for ASE.
    template_path = 'ZnO_template.in'  
   
   # (1,1,1) to (3,3,3) supercells
    #for nx in range(1, 4):  
    #    for ny in range(1, 4):  
    #        for nz in range(1, 4):
    #            make_pw_input_from_existing_in(input_file_path=template_path,
    #                a=CELLDM1_Angstroms,
    #                covera=CELLDM3,
    #                cellsize=(nx, ny, nz) )
    
    # ================ STRAIN =================== #
    for a in a_values:
        for covera in c_values:
            make_pw_input_from_existing_in(input_file_path=template_path,
                a=a,
                covera=covera,
                cellsize=(1, 1, 3))

