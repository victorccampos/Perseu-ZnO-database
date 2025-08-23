#!/usr/bin/env python3
import subprocess
from pathlib import Path
from datetime import datetime
from itertools import product

QE_BINARY = "pw.x"
INPUT_FILE = Path("supercells/ZnO-111.in")
OUTPUT_DIR = Path("benchmarks")
OUTPUT_DIR.mkdir(exist_ok=True)


class BenchmarkConfig:
    def __init__(self, nproc=None, nk_values=None, nb_values=None, nt_values=None, nd_values=None, mode="flags"):
        """
        mode = "procs" | "flags" | "both"
        """
        self.nproc = nproc
        self.nk_values = nk_values
        self.nb_values = nb_values
        self.nt_values = nt_values
        self.nd_values = nd_values
        self.mode = mode


def run_qe_job(input_path: Path, output_path: Path,
               num_processes: int, nk: int, nb: int, nt: int, nd: int) -> None:
    """
    Executa um único job do QE com paralelização controlada.
    """
    command_line_text = (
        f"mpirun -np {num_processes} {QE_BINARY} "
        f"-nk {nk} -nb {nb} -nt {nt} -nd {nd} "
        f"-in {str(input_path)}"
    )
    command_line = command_line_text.split()

    with output_path.open(mode="w") as f_out:
        subprocess.run(command_line, check=True, stdout=f_out, stderr=subprocess.STDOUT)

    print(f"✅ Job finalizado: {output_path.name} -- {datetime.now().strftime('%H:%M:%S')}")


def generate_scenarios(config: BenchmarkConfig):
    """
    Gera  as combinações de benchmark de acordo com o modo.
    """
    if config.mode == "procs":
        for nproc in config.nproc:
       (nproc, 1, 1, 1, 1)  # apenas varia nº de processos

    elif config.mode == "flags":
        nproc = config.nproc[0]
        for nk, nb, nt, nd in product(config.nk_values, config.nb_values,
                                      config.nt_values, config.nd_values):
            if nk * nb * nt * nd == nproc:
  (nproc, nk, nb, nt, nd)

    elif config.mode == "both":
        for nproc in config.nproc:
            for nk, nb, nt, nd in product(config.nk_values, config.nb_values,
                                          config.nt_values, config.nd_values):
                if nk * nb * nt * nd == nproc:
                    yield (nproc, nk, nb, nt, nd)


def run_benchmark(config: BenchmarkConfig):
    """
    Executa todos os cenários de benchmark.
    """
    for nproc, nk, nb, nt, nd in generate_scenarios(config):
        out_file = OUTPUT_DIR / f"bench_p{nproc}_nk{nk}_nb{nb}_nt{nt}_nd{nd}.out"
        try:
            run_qe_job(INPUT_FILE, out_file, nproc, nk, nb, nt, nd)
        except subprocess.CalledProcessError:
            print(f"❌ Falhou: procs={nproc}, nk={nk}, nb={nb}, nt={nt}, nd={nd}")


if __name == "__main__":
    # Exemplo 1: variar só flags, com 64 processos fixos
    config = BenchmarkConfig(
        nproc=[64],
        nk_values=[1, 2, 4, 8],
        nb_values=[1, 2, 4, 8],
        nt_values=[1, 2, 4, 8],
        nd_values=[1, 2, 4, 8],
        mode="flags"
    )

    # Exemplo 2: variar só número de processos
    # config = BenchmarkConfig([16, 32, 64], [1], [1], [1], [1], mode="procs")

    # Exemplo 3: variar processos e flags juntos
    # config = BenchmarkConfig([32, 64], [1, 2, 4], [1, 2], [1, 2], [1, 2], mode="both")

    run_benchmark(config)
