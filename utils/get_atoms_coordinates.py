import argparse
import re
import pandas as pd
from pathlib import Path
import pprint


parser = argparse.ArgumentParser(description="obtém coordenadas atômicas de pw.scf.in")
parser.add_argument("-f", "--file", help="Arquivo de input para se obter as coordenadas")

args = parser.parse_args()

filepath_arg: str = args.file

def ler_posicoes(filepath: str) -> pd.DataFrame:
    """ Obtém DataFrame com posições iniciais do input pw scf"""
    with open(filepath, 'r') as file:
        texto = file.read()

        # Encontra a seção ATOMIC_POSITIONS
        match = re.search(r'ATOMIC_POSITIONS.*?\n(.*?)(\n\n|\Z)', texto, re.S)
        if not match:
            raise ValueError(f'Secção ATOMIC_POSITIONS não encontrada em {filepath}')

        bloco = match.group(1).strip().splitlines()

        # Parseia linhas em dataframe
        dados = []
        for linha in bloco:
            # print(linha)
            partes = linha.split()
            elemento = partes[0]
            x, y, z = map(float, partes[1:4])
            
            dados.append({'elemento': elemento, 'x': x, 'y': y, 'z': z})    
        
        
        return pd.DataFrame(dados)


if __name__ == '__main__':
    df_posicoes_n_pertubadas = ler_posicoes(filepath=filepath_arg)

    # Processar todos outputs.
    input_dir = Path('./SCFs_rnd.in')

    resultados_variacoes = []

    for input_file in sorted(input_dir.glob('scf-*.in')):
        df_in = ler_posicoes(input_file)


        # Cálculo da variação percentual |delta| / valor do template
        variacao = (df_in[['x', 'y', 'z']] - df_posicoes_n_pertubadas[['x', 'y', 'z']]).abs() / df_posicoes_n_pertubadas[['x', 'y', 'z']] * 100

        variacao['arquivo'] = input_file.name
        variacao['elemento'] = df_in['elemento']

        resultados_variacoes.append(variacao)

    pprint.pprint(resultados_variacoes[:5])

