# %%
import re
from pathlib import Path
import pandas as pd
import numpy as np


def ler_posicoes(filepath):
    """Lê as posições atômicas a partir de um arquivo .in ou template."""
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
        print(linha)
        partes = linha.split()
        elemento = partes[0]
        x, y, z = map(float, partes[1:4])
        
        dados.append({'elemento': elemento, 'x': x, 'y': y, 'z': z})

    return pd.DataFrame(dados)


# 🔧 Configurações
template_path = Path('./template.in')
input_dir = Path('./SCFs_rnd.in')  # ajuste conforme seu diretório


# 📥 Ler template
df_template = ler_posicoes(template_path).reset_index(drop=True)

# %%
# 📦 Processar todos os arquivos .in
resultados = []

for input_file in sorted(input_dir.glob('*.in')):
    df_in = ler_posicoes(input_file).reset_index(drop=True)

    if len(df_template) != len(df_in):
        raise ValueError(f'Número de átomos diferente em {input_file}')

    # Cálculo da variação percentual |delta| / valor do template
    variacao = (df_in[['x', 'y', 'z']] - df_template[['x', 'y', 'z']]).abs() / df_template[['x', 'y', 'z']] * 100

    variacao['arquivo'] = input_file.name
    variacao['elemento'] = df_in['elemento']

    resultados.append(variacao)


# 🔗 Concatenar todos os resultados
df_variacoes = pd.concat(resultados, ignore_index=True)


# 📊 Estatísticas gerais
estatisticas = df_variacoes[['x', 'y', 'z']].describe()

print('==== Variações Percentuais nas Coordenadas ====')
print(estatisticas)

# 🔥 Variação máxima observada
max_var = df_variacoes[['x', 'y', 'z']].max().max()
print(f'\n⚠️ Variação percentual máxima encontrada: {max_var:.6f}%')

# ✔️ Salvar em CSV, se desejar
#df_variacoes.to_csv('variacoes_atomicas.csv', index=False)
