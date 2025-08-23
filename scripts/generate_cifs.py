from pathlib import Path
from ase.io import read
from typing import Optional

def get_cif(input_file: Path, output_dir: Optional[Path] = None) -> None:
    """Gera arquivo '.cif' a partir arquivo '.in'"""
    
    if not input_file.is_file():
        print(f"Erro: O arquivo de entrada '{input_file}' não foi encontrado.")
        return
    
    cif_name = input_file.name.replace(".in", ".cif")
    
    output_path = output_dir / cif_name if output_dir else Path(cif_name)
    
    try:
        atoms = read(input_file, format="espresso-in")
        atoms.write(output_path, format="cif") 
        
        print(f"Arquivo {cif_name} gerado.")
    except Exception as e:
        print(f"Falha ao processar o arquivo '{input_file.name}': {e}")


if __name__ == "__main__":
    modo_diretorio = True
    # MODO DIRETORIO
    diretorio_dos_inputs = Path('./cell323.in') # TODO: Mudar diretórios de forma mais fácil.
    diretorio_de_saida =  Path('./cifs_' + diretorio_dos_inputs.stem)

    # MODO INPUT ÚNICO
    single_input_path = Path(f'./cell323.in/ZnO-3.40-1.58-323.in')

    if modo_diretorio:
        diretorio_de_saida.mkdir(parents=True, exist_ok=True)
        print("Gerando CIF - DIRÉTORIO")
        print(f'\t {diretorio_dos_inputs}')
        input_files = sorted([file for file in diretorio_dos_inputs.iterdir() if file.is_file() and file.suffix == ".in"])

        for file in input_files:
            get_cif(input_file=file, output_dir=diretorio_de_saida)

        print(f"\nConversão finalizada. Arquivos salvos em '{diretorio_de_saida}'.")
        
    else:
        get_cif(single_input_path)
    
    
