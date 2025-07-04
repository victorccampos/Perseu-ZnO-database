from ase.build import bulk
from ase.io.espresso import write_espresso_in


# Dicionário com seções
# &control, &system, &electrons

header_of_input = {
    'control': {'calculation': 'scf',
                'prefix': 'ZnO_supercells',
                'pseudor_dir': r'../../../pseudos',
                'outdir': './', 'tprnfor': True },
    'system': {'ibrav': 0, 'ecutwfc': 50, 'ecutrho': 400, 'occupations': 'fixed'},
    'electrons': {'conv_thr':1.0e-8, 'mixing_beta': 0.3}
}

pseudopotentials_dict = {
    'Zn': 'Zn.upf',
    'O': 'O.upf'
}


# Criar arquivo de input.
write_espresso_in(file='ZnO-222.in',
                  atoms = supercell,
                  input_data = input_header,
                  kpts=(6,6,6),
                  pseudopotentials=pseudo_potentials
)
