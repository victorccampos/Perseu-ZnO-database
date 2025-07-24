from pathlib import Path
import subprocess
import argparse
import sys
import os


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments for the QE input runner."""
    
    parser = argparse.ArgumentParser(
        description="Run all Quantum ESPRESSO input files in a specified directory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-d",
        "--directory",
        type=str,
        default="SFCs.in",
        help="Directory containing Quantum ESPRESSO input files (.in)."
    )
    parser.add_argument(
        "-np",
        "--num-processes",
        type=str,
        default=16,
        help="Number of MPI processes to use for each QE run."
    )
    parser.add_argument("-nk",
        "--npools",
        type=str,
        default=2,
        help="Number of pools for FFT parallelization."
    )
    return parser.parse_args()


def run_qe_job(input_path: Path, output_path: Path, num_processes: str, npools: str) -> None:
    """
    Execute a single Quantum ESPRESSO job for the given input file.
    Args:
        num_processes (int): Number of MPI processes to use.
        npools (int): Number of pools for FFT parallelization.
    """
    
    cmd_str = f"mpirun -np {num_processes} pw.x -npools {npools} -in {str(input_path)}"
    cmd: list[str] = cmd_str.split()

    try:
        with output_path.open(mode='w') as f_out:
            subprocess.run(
                cmd,
                check=True,
                stdout=f_out,
                stderr=subprocess.STDOUT,
            )
        print(f"Successfully completed: {input_path.name}")

    except subprocess.CalledProcessError as e:
        print(f"Failed to execute {input_path.name}: {e}")
    except subprocess.TimeoutExpired:
        print(f"Timeout expired for {input_path.name}")
    except Exception as e:
        print(f"Unexpected error for {input_path.name}: {e}")

def main():
    
    args = parse_arguments()
    
    # If hybrid parallelization is used (OpenMP + MPI)
    os.environ['OMP_NUM_THREADS'] = '4'    
    
    # Setup directories
    base_dir = Path(__file__).parent
    input_dir = base_dir / args.directory
    output_dir = base_dir / args.directory.replace(".in", ".out")
    
    try:
        output_dir.mkdir(exist_ok=True)
    except Exception as e:
        print(f"Failed to create output directory {output_dir}: {e}")
        sys.exit(1)

    
    try:
        input_files = sorted(
            [file for file in input_dir.iterdir()
             if file.is_file() and file.name.endswith(".in")]
        )
        if not input_files:
            print(f"No .in files found in {input_dir}")
            sys.exit(1)
    except Exception as e:
        print(f"Error accessing input directory {input_dir}: {e}")
        sys.exit(1)
    
    # Process files sequentially
    for n, input_path in enumerate(input_files):
        output_filename = input_path.stem + ".out"
        output_path = output_dir / output_filename
        run_qe_job(input_path, output_path, args.num_processes, args.npools)

    
if __name__ == "__main__":
    main()
