import sys
import re


# Output file from self-consistent field calculation
output_file = sys.argv[1]

with open(output_file, 'r') as f:
    data = f.read()

#========================== Extract total energy (Ry) ========================#
TOTAL_ENERGY_PATTER_IN_QE = r"!\s+total energy\s+=\s+([-\d.]+)\s+Ry"
energy_match = re.search(
    pattern = TOTAL_ENERGY_PATTER_IN_QE,
    string = data)
total_energy = float(energy_match.group(1)) if energy_match else 0.0

#========================= Extract forces in (Ry/au) ==========================#
                                    # Forces acting on atoms (cartesian axes, Ry/au):
# atom    1 type  1   force =       # non-local contrib.
# atom    2 type  1   force =       # ionic contribution
# atom    3 type  2   force =       # local contribution
# atom    4 type  2   force =       # core correction contribution
                                    # Hubbard contrib.
                                    # SCF correction ter

FORCES_PATTERN_in_QE = r"atom\s+\d+\s+type\s+\d+\s+force\s+=\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"

forces = []
atomic_species = []

force_matches = re.finditer(FORCES_PATTERN_in_QE, data)
i_trunc = 4 # Primeiro bloco de forças têm 4 linhas no PWSCF do QE

for index, match in enumerate(force_matches):
    if index < i_trunc:
	    forces.append([
        float(match.group(1)), # Fx
        float(match.group(2)), # Fy
        float(match.group(3))  # Fz
        ])
	    atomic_species.append(match.group().split[:4])
# pprint.pprint(forces)
# =============================================================================#



# Escreve no formato RuNNer (usado pelo Ænet)
with open("train.data", "a") as f:
    f.write(f"begin\n")
    f.write(f"comment ZnO structure\n")
    f.write(f"lattice __CELLDM1__ 0.0 0.0\n")
    f.write(f"lattice 0.0 __CELLDM1__ 0.0\n")
    f.write(f"lattice 0.0 0.0 __CELLDM3__\n")
    f.write(f"energy {total_energy}\n")
    for i, force in enumerate(forces):
        f.write(f"atom {force[0]} {force[1]} {force[2]} 0.0 0.0 0.0 Zn\n")  # Ajuste para O se necessário
    f.write(f"end\n")
