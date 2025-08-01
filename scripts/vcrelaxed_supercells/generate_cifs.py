import glob
from ase.io import read


# Arquivos de input scf
input_files = glob.glob("ZnO-*.in")

for file in input_files:
    atoms = read(file, format="espresso-in")
    cif_name: str = file.replace("ZnO-", "cell-").replace(".in", ".cif")

    # Salvando em formato cif com ASE.
    atoms.write(cif_name, format="cif") # format = ["vasp", "xyz", "json"]
    print(f"Arquivo {cif_name} gerado.")