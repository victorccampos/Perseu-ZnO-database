import numpy as np
import os
from typing import Tuple

# Atomic Simulation Environment (ASE) imports
from ase.io import read                         # read ZnO_template.in
from ase import Atoms
from ase.build import make_supercell
from ase.io.espresso import write_espresso_in



def make_pw_from_template(
    template_path: str,
    a: float,
    covera: float,
    cellsize: Tuple[int, int, int]
    ) -> None:
    """Create a supercell input file for Quantum ESPRESSO based on an existing input file."""
    
    
    primitive_cell = read(template_path, format='espresso-in') 
    nx, ny, nz = cellsize
    supercell_matrix = np.array([
        [nx, 0, 0],
        [0, ny, 0],
        [0, 0, nz] 
    ])
    supercell = make_supercell(prim=primitive_cell, P=supercell_matrix)
    
    
    cell_dirname = f'cell{nx}{ny}{nz}.in'
    os.makedirs(cell_dirname, exist_ok=True)

    # Input file name based on the parameters
    input_name: str = f'ZnO-{a:.2f}-{covera:.2f}-{nx}{ny}{nz}.in'
    
    print(f'Produzindo : {input_name}')
    
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
    KPOITNS_ORIGINAL = np.array([6, 6, 6])  
    KPOINT_SUPERCELL = KPOITNS_ORIGINAL // np.array(cellsize)  
    KPOINT_SUPERCELL[KPOINT_SUPERCELL == 0] = 1 

    
    write_espresso_in(
        file = os.path.join(cell_dirname, input_name), 
        atoms = supercell,
        input_data = HEADER_INPUT,
        pseudopotentials = PSEUDOS,
        kpts = tuple(KPOINT_SUPERCELL),
        koffset = (0,0,0),
        crystal_coordinates = True # TODO: check if this is needed
    )


if  __name__ == "__main__":
    percent_range = np.arange(start= -0.10, stop=0.12, step=0.02)

    # Relaxed structure parameters for ZnO
    CELLDM1_Bohr_units = 6.178821408099141            
    CELLDM1_angstroms = CELLDM1 * 0.52918  
    CELLDM3 = 1.614358356153010 

    # Strained parameters - FUTURE
    a_values = CELLDM1_Angstroms * (1 + percent_range)
    c_values = CELLDM3 * (1 + percent_range)
    
    # Path to the template input file - 1x1x1 atoms structure for ASE.
    template_path = 'ZnO_template.in'  
   
    
    # ================ STRAIN =================== #
    for a in a_values:
        for covera in c_values:
            make_pw_from_template(
                template_path=template_path,
                a=a,
                covera=covera,
                cellsize=(2, 2, 2)
                )

