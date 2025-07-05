import numpy as np
from ase.build import make_supercell, bulk
from ase.io.espresso import write_espresso_in
from typing import Tuple



#======================================================================================#

percent_range = np.arange(start= -0.10, stop=0.12, step=0.02)
# [0.9  0.92 0.94 0.96 0.98 1.   1.02 1.04 1.06 1.08 1.1 ]

CELLDM1 = 6.178821408099141 # in Bohr units
CELLDM3 = 1.614358356153010 # Adimensional c/a ratio

CELLDM1_Angstroms = CELLDM1 * 0.52918

a_values = CELLDM1_Angstroms * percent_range

# Input comum para todas as simulações

#======================================================================================#


def make_pw_input(a: float, covera: float, cellsize: Tuple[int, int, int]):
    primitive_cell = bulk(name='ZnO', crystalstructure='wurtzite', a=a, covera=covera)
    
    nx, ny, nz = cellsize
    supercell_matrix = np.array([
        [nx, 0, 0],
        [0, ny, 0],
        [0, 0, nz] ]
    )
    supercell = make_supercell(prim=primitive_cell, P=supercell_matrix)
    
    
    
    input_name = f'ZnO-{a:.2f}-{covera:.2f}-{nx}{ny}{nz}.in'
    print(f'Produzindo : {input_name}')
    
    HEADER_INPUT = {
    'control': {
        'calculation': 'scf',
        'prefix': f'{input_name}',
        'pseudo_dir': '../pseudos',
        'outdir': './',
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
    
    # Pseudopotenciais
    PSEUDOS = {
        'Zn': 'Zn.upf',
        'O': 'O.upf'
    }
    KPOINTS = (6,6,6)
    
    write_espresso_in(
        file = input_name,
        atoms = supercell,
        input_data = HEADER_INPUT,
        pseudopotentials = PSEUDOS,
        kpts = KPOINTS,
        koffset = (0,0,0),
        crystal_coordinates = True
    )


make_pw_input(a=3.269691715, covera=CELLDM3, cellsize=(1, 1, 1))


