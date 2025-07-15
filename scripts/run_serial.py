from pathlib import Path
import subprocess
import argparse


# ================================= PARSER ==================================== #

parser = argparse.ArgumentParser(
        description="Rodar todos inputs QE contidos em um diretório.")
    
# Folder Flag
parser.add_argument("-d", "--directory",
                    default="SFCs.in",
                    help="Nome diretório onde estão armazenados os inputs do QE.")

args = parser.parse_args()
# Pressupõe que os inputs já foram gerados.
# ============================================================================= #

# Diretórios 
BASE_DIR = Path(__file__).parent

input_dir = BASE_DIR / args.directory
output_dir = BASE_DIR / args.directory.replace(".in", ".out")
output_dir.mkdir(exist_ok=True)     # Cria se não existir


input_files = sorted(
    [file for file in input_dir.iterdir()
     if file.is_file() 
     and file.name.endswith(".in")]
) 


for n, input_path in enumerate(input_files):
    print(f'Executando script {n+1}:'+'\n\t'+f'{input_path.name}')
    
    
    output_filename = input_path.stem + ".out" # scf.{}.{} + .out
    output_path = output_dir / output_filename
    print('Saída será salva em:\n\t', output_path)
    

    # Abre o arquivo de saída e redireciona stdout e stderr para ele
    with output_path.open(mode='w') as f_out:
        try:
            subprocess.run(
                # Executa o comando mpirun
                ["mpirun", "-np", "8", "pw.x", "-in", str(input_path)],
                check=True,
                stdout=f_out ,
                stderr=subprocess.STDOUT
            )
        except subprocess.CalledProcessError as e:
            print(f"Erro ao executar {input_path.name}: {e}")
            continue # tenta o próximo input se erro
    
