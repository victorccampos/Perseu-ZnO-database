import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
from pathlib import Path

# Diretórios com outputs
diretorio_de_outputs: str = '../SCFs_rnd.out/' 

output_dir = Path(diretorio_de_outputs)

def get_energy_dataframe(output_dir: Path, ) -> pd.DataFrame:
    arquivos = sorted(
        output_dir.glob('*.out')
    )

    dados = []

    for arquivo in arquivos:
        pass
    
    return []


