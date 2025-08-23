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

    primitive_cell = read(template_path, format='espresso-in')
    nx, ny, nz = supercell_size
    supercell_matrix = np.diag(supercell_size)
    supercell = make_supercell(prim=primitive_cell, P=supercell_matrix)

    # Base names
    cell_dirname = f'cell{nx}{ny}{nz}.in'
    input_name: str = f'ZnO-{a:.2f}-{covera:.2f}-{nx}{ny}{nz}.in'

    # Random Displacements : rattle()
    is_noisy: bool = add_noise and noise_stdev > 0
    seed = np.random.randint(0, 1e9) 

    if is_noisy:
        cell_dirname = f'cell{nx}{ny}{nz}-{noise_stdev}.in'
        input_name = f'ZnO-{a:.2f}-{covera:.2f}-{nx}{ny}{nz}-{noise_stdev}.in'
        supercell.rattle(stdev=noise_stdev, seed=seed)


    print(f'Produzindo input self-consistent-field == {input_name}')
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

    # Adapt K-grid to larger supercells, mantaining grid density
    KPOITNS_ORIGINAL = np.array([6, 6, 6])
    KPOINT_SUPERCELL = KPOITNS_ORIGINAL // np.array(supercell_size)
    KPOINT_SUPERCELL[KPOINT_SUPERCELL == 0] = 1


    write_espresso_in(
        file = os.path.join(cell_dirname, input_name),
        atoms = supercell,
        input_data = HEADER_INPUT,
        pseudopotentials = PSEUDOS,
        kpts = tuple(KPOINT_SUPERCELL),
        koffset = (0,0,0),
        crystal_coordinates = True
    )

if  __name__ == "__main__":
    template_path = 'ZnO_template.in'
    
    percent_range = np.arange(start= -0.10, stop=0.12, step=0.02)

    # Relaxed structure parameters for ZnO
    CELLDM1_Bohr_units = 6.178_821_408_099_141
    CELLDM1_angstroms = CELLDM1_Bohr_units * 0.52918
    CELLDM3 = 1.614_358_356_153_010


    strained_a_values = CELLDM1_angstroms * (1 + percent_range)
    strained_covera_values = CELLDM3 * (1 + percent_range)
    
    # ================ STRAIN =================== #
    noise_level = 0.05

    for a in strained_a_values:
        for covera in strained_covera_values:
            make_scf_from_template(template_path, a, covera,(1, 1, 2), add_noise=True, noise_stdev=noise_level)

