import numpy as np
import os
from typing import Tuple

# Atomic Simulation Environment (ASE) imports
from ase.io import read                         
from ase import Atoms
from ase.build import make_supercell
from ase.io.espresso import write_espresso_in


def make_scf_from_template(
    template_path: str,
    a: float,
    covera: float,
    supercell_size: Tuple[int, int, int],
    add_noise: bool = False,
    noise_stdev: float = 0.015
    ) -> None:
    """Create a supercell input file for Quantum ESPRESSO based on an existing input file.

    Args:
        a : lattice paremeter in Angstrom
        covera: ratio c/a
        noise_stdev: Desloca os átomos com uma distribuição gaussiana cujo desvio padrão é de noise_stdev Å

        For more info about random displacement: https://ase-lib.org/_modules/ase/atoms.html#Atoms.rattle
    """
    # Lê a estrutura do template para obter a base de átomos (1x1x1)
    primitive_cell_template = read(template_path, format='espresso-in')
    
    ############################################################################
    # CELL_PARAMETERS (v1,v2,v3)
    c = a * covera
    new_cell_vectors = [
        [a, 0, 0],
        [-a/2.0, a*(np.sqrt(3)/2.0), 0 ], 
        [0, 0, c]
    ]
    
    new_cell: Atoms = Atoms(
        symbols = primitive_cell_template.get_chemical_symbols(),
        scaled_positions = primitive_cell_template.get_scaled_positions(),
        cell = new_cell_vectors,
        pbc = True
    )
    ############################################################################
    
    nx, ny, nz = supercell_size
    supercell_matrix = np.diag(supercell_size)
    
    # Criando supercélula a partir da nova célula 
    supercell = make_supercell(prim=new_cell, P=supercell_matrix)

    # Base names
    cell_dirname = f'cell{nx}{ny}{nz}.in'
    input_name: str = f'ZnO-{a:.2f}-{covera:.2f}-{nx}{ny}{nz}.in'

    # Random Displacements : rattle()
    is_noisy: bool = add_noise and noise_stdev > 0
    seed = np.random.randint(0, 1e9) 

    if is_noisy:
        cell_dirname = f'random_cell{nx}{ny}{nz}-{noise_stdev}.in'
        input_name = f'ZnO-{a:.2f}-{covera:.2f}-{nx}{ny}{nz}-{noise_stdev}.in'
        supercell.rattle(stdev=noise_stdev, seed=seed)
    # TODO: Quando adicionar ruído, para estrutura fixa, colocar número da estrutura. 
    # Ex: ZnO-{a:.2f}-{covera:.2f}-{nx}{ny}{nz}-{noise_stdev}-{i}.in

    print(f'Produzindo input self-consistent-field => {input_name}')
    os.makedirs(cell_dirname, exist_ok=True)

    # Input SCF Quantum Espresso.
    prefix_qe = input_name.replace(".in", "")
    HEADER_INPUT = {
    'control': {
        'calculation': 'scf',
        'prefix': prefix_qe,
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
    K_GRID = (6, 6, 6) 
    
    # Generate ZnO-{a}-{c/a}.in or ZnO-{a}-{c/a}-{noise_level}.in 
    write_espresso_in(
        file = os.path.join(cell_dirname, input_name),
        atoms = supercell,
        input_data = HEADER_INPUT,
        pseudopotentials = PSEUDOS,
        kpts = K_GRID,
        koffset = (0,0,0),
        crystal_coordinates = True
    )

if  __name__ == "__main__":
    template_path = 'ZnO_template.in'
    
    percent_range = np.arange(start= -0.10, stop=0.12, step=0.02) # [start,stop)
    print(f"Valores de stain (%): {percent_range}")
    # Relaxed structure parameters for ZnO
    # lattice parameter a
    celldm1_bohr = 6.178_821_408_099_141
    celldm1_angstroms = celldm1_bohr * 0.52918
    
    # ratio c/a
    celldm3 = 1.614_358_356_153_010

    strained_a_values = celldm1_angstroms * (1 + percent_range)
    strained_covera_values = celldm3 * (1 + percent_range)
    
    # noise_level = 0.05
    # ================ STRAIN LOOP =================== #
    
    for a in strained_a_values:
        for covera in strained_covera_values:
            make_scf_from_template(template_path, a, covera, supercell_size=(3, 3, 1))

