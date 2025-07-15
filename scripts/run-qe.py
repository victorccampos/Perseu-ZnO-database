from pathlib import Path
import subprocess
import argparse
import logging
from typing import List
import sys
import os

# ============================== CONFIG LOGGER ================================ #

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('qe_runner.log')
    ]
)
logger = logging.getLogger(__name__)

# ============================== PARSER SETUP ================================ #

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for the QE input runner.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
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
        type=int,
        default=16,
        help="Number of MPI processes to use for each QE run."
    )
    parser.add_argument(
        "--npool",
        type=int,
        default=2,
        help="Number of pools for FFT parallelization."
    )
    return parser.parse_args()

# ============================ MAIN PROCESSING ================================ #

def run_qe_job(input_path: Path, output_path: Path, num_processes: int, npool: int) -> None:
    """
    Execute a single Quantum ESPRESSO job for the given input file.

    Args:
        input_path (Path): Path to the QE input file.
        output_path (Path): Path where the output file will be saved.
        num_processes (int): Number of MPI processes to use.
        npool (int): Number of pools for FFT parallelization.
        pw_path (str): Path to the pw.x executable.
    """
    logger.info(f"Executing QE job: {input_path.name} (Output: {output_path})")
    cmd = [
        "mpirun", "-np", str(num_processes), "--map-by", "numa:pe=2",
        "pw.x", "-in", str(input_path), "-npool", str(npool)
    ]
    try:
        with output_path.open(mode='w') as f_out:
            subprocess.run(
                cmd,
                check=True,
                stdout=f_out,
                stderr=subprocess.STDOUT,
            )
        logger.info(f"Successfully completed: {input_path.name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to execute {input_path.name}: {e}")
    except subprocess.TimeoutExpired:
        logger.error(f"Timeout expired for {input_path.name}")
    except Exception as e:
        logger.error(f"Unexpected error for {input_path.name}: {e}")

def main() -> None:
    """
    Main function to orchestrate the execution of Quantum ESPRESSO jobs.

    Sets up directories, collects input files, and processes them sequentially.
    """
    args = parse_arguments()
    
    # Set OMP_NUM_THREADS environment variable
    os.environ["OMP_NUM_THREADS"] = "2"
    
    # Setup directories
    base_dir = Path(__file__).parent
    input_dir = base_dir / args.directory
    output_dir = base_dir / args.directory.replace(".in", ".out")
    
    try:
        output_dir.mkdir(exist_ok=True)
    except Exception as e:
        logger.error(f"Failed to create output directory {output_dir}: {e}")
        sys.exit(1)

    # Collect input files'
    try:
        input_files = sorted(
            [file for file in input_dir.iterdir()
             if file.is_file() and file.name.endswith(".in")]
        )
        if not input_files:
            logger.error(f"No .in files found in {input_dir}")
            sys.exit(1)
    except Exception as e:
        logger.error(f"Error accessing input directory {input_dir}: {e}")
        sys.exit(1)

    logger.info(f"Found {len(input_files)} input files to process.")
    
    # Process files sequentially
    for n, input_path in enumerate(input_files):
        logger.info(f'Executing job {n+1}/{len(input_files)}: {input_path.name}')
        output_filename = input_path.stem + ".out"
        output_path = output_dir / output_filename
        logger.info(f'Output will be saved to: {output_path}')
        run_qe_job(input_path, output_path, args.num_processes, args.npool)

    logger.info("All QE jobs completed.")

if __name__ == "__main__":
    main()
