from pathlib import Path
import subprocess
import os
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed


# ============================ PARSER ====================================== #

parser = argparse.ArgumentParser(
    description="Rodar todos inputs QE contidos em um diretório em paralelo."
)

parser.add_argument(
    "-d", "--directory",
    default="SCFs.in",
    help="Nome do diretório onde estão armazenados os inputs do QE."
)

parser.add_argument(
    "-np", "--nprocess",
    default=8,
    type=int,
    help="Número de processadores por cálculo (padrão=8)"
)

args = parser.parse_args()

# ============================ CONFIGURAÇÕES =============================== #

# Diretórios
BASE_DIR = Path(__file__).parent

input_dir = BASE_DIR / args.directory
output_dir = BASE_DIR / args.directory.replace(".in", ".out")
output_dir.mkdir(exist_ok=True)

# Arquivos de input
input_files = sorted([
    file for file in input_dir.iterdir()
    if file.is_file() and file.name.startswith("scf-") and file.name.endswith(".in")
])


TOTAL_CPUS = os.cpu_count() # 64
CPUS_PER_JOB = args.nprocess             # default 8
MAX_WORKERS = TOTAL_CPUS // CPUS_PER_JOB # default 8


# ============================ FUNÇÃO DE EXECUÇÃO ========================== #

def run_job(input_path):
    output_filename = input_path.stem + ".out"
    output_path = output_dir / output_filename

    print(f"\nExecutando {input_path.name} com saída em {output_filename}")

    with output_path.open(mode='w') as f_out:
        try:
            subprocess.run(
                ["mpirun", "-np", str(CPUS_PER_JOB), "../bin/pw.x", "-in", str(input_path)],
                check=True,
                stdout=f_out,
                stderr=subprocess.STDOUT
            )
            print(f"Finalizado {input_path.name}")
        except subprocess.CalledProcessError as e:
            print(f"Erro ao executar {input_path.name}: {e}")


# ============================ EXECUTOR ==================================== #

if __name__ == "__main__":
    print(f"\nRodando com até {MAX_WORKERS} jobs simultâneos, cada um usando {CPUS_PER_JOB} CPUs.\n")

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [executor.submit(run_job, input_file) for input_file in input_files]

        for future in as_completed(futures):
            future.result()  # Levanta exceções, se houver

    print("\n=== TODOS OS JOBS FINALIZADOS ===")
