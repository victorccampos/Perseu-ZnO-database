from ase.build import bulk
from ase.io.espresso import write_espresso_in
from typing import Tuple

# Parâmetros da estrutura primitiva de ZnO (wurtzita)
a = 3.25  # parâmetro de rede a (angstrom)
c = 5.2   # parâmetro de rede c (angstrom)

# Pseudopotenciais para Zn e O (nomes devem coincidir com os arquivos no pseudo_dir)
pseudopotentials = {
    'Zn': 'Zn.upf',
    'O': 'O.upf'
}

# Input comum para todas as simulações
input_data = {
    'control': {
        'calculation': 'scf',
        'prefix': 'ZnO_supercell',
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

# Função para gerar supercélula e input QE
def gerar_input_qe(cellsize: Tuple[int, int, int], a, c_over_a):
    """
    Params:
    cellsize: tuple ex: 1,1,1; 2,1,1
    a: lattice parameter
    c_over_a: lattice parameter to ase.bulk() function
    """
    nx, ny, nz = cellsize
    label = f"{nx}{ny}{nz}"
    
    # Célula primitiva e supercélula
    atoms = bulk('ZnO', 'wurtzite', a=a, covera=c_over_a)
    supercell = atoms.repeat((nx, ny, nz))

    # Prefixo e nome do arquivo com base na supercélula
    input_data['control']['prefix'] = f'ZnO_{label}'
    filename = f'ZnO-{label}.in'

    # Exporta input QE com K-points fixos
    write_espresso_in(
        file=filename,
        atoms=supercell,
        input_data=input_data,
        pseudopotentials=pseudopotentials,
        kpts=(6, 6, 6)
    )
    print(f"[✓] Gerado: {filename}")

# Loop para várias supercélulas
celldm1_in_angstrom = 3.269_708_713
celldm3_in_angstrom = 1.614_358_356

for nx in range(1, 4):
    for ny in range(1, 4):
        for nz in range(1, 4):
            gerar_input_qe(cellsize=(nx, ny, nz), a=celldm1_in_angstrom, c_over_a=celldm3_in_angstrom)
