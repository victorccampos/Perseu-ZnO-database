from pathlib import Path
import subprocess
import argparse
import os
from datetime import datetime

# (8, 2) -> 16 cores
# (16, 2) -> 32 cores
# (24, 2) -> 48 cores
# (32, 2) -> 64 cores
# (40, 1) -> 40 cores
# (48, 1) -> 48 cores
# (56, 1) -> 56 cores
# (64, 1) -> 64 cores

# Arquivos e Diretórios
# BASE_DIR = Path(__file__).parent
# cellsize = '333'
# input_file: str = f"ZnO-{cellsize}.in" 

# # CPU - os.environ['OMP_NUM_THREADS']
# total_cores = 64

# procs_mpi = [nproc for nproc in range(4, total_cores+1, 4)]

# for i in procs_mpi:
#     output_name = f"bench_{cellsize}_{i}.out"
#     #nthreads = max(1, min(2, total_cores // i))
#     nthreads = '1'
    
#     #os.environ['OMP_NUM_THREADS'] = '1'
#     #os.environ['OMP_NUM_THREADS'] = str(nthreads)
#     print(f"Executando {output_name}: {i} processos MPI x {nthreads} threads OMP = {i * nthreads} cores")
    

#     command_line_text = f"mpirun -np {i} pw.x -in {input_file}"    
#     command_line: list[str] = command_line_text.split()
#     output_path = BASE_DIR / output_name
    
    
#     with output_path.open(mode='w') as f_out:
#         try:
#             subprocess.run(
#                 command_line,
#                 check=True,
#                 stdout=f_out ,
#                 stderr=subprocess.STDOUT
#             )
#         except subprocess.CalledProcessError as e:
#             print(f"Erro ao executar {input_file}: {e}")
#             continue 


def run_qe_job(input_path: Path, output_path: Path, num_processes: str, npools: str) -> None:
    """
    Execute a single Quantum ESPRESSO job for the given input file.
    Args:
        num_processes (int): Number of MPI processes to use.
        npools (int): Number of pools for FFT parallelization.
    """

    command_line_text = f"mpirun -np {num_processes} pw.x -npools {npools} -in {str(input_path)}"
    command_line: list[str] = command_line_text.split()
    
    # input_path não deve conter espaços!
    with output_path.open(mode='w') as f_out:
            subprocess.run(command_line, check=True, stdout=f_out, stderr=subprocess.STDOUT)
        # Log nohup.out
    print(f"Job completed: {input_path.name} -- {datetime.now().strftime('%H:%M:%S')}")




if __name__ == "__main__":
     supercell_in: str = "ZnO-333.in"
     input_path = Path(__file__).parent / supercell_in 
     output_filename = input_path.stem + f'{}.out'   
     print()
  
     
   
    
