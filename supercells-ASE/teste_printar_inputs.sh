#!/bin/bash

# Diretório atual
INPUT_DIR="."

# Obtém todos os arquivos .in no diretório atual, ordenados numericamente
INPUT_FILES=($(ls -v "$INPUT_DIR"/*.in 2>/dev/null))

# Verifica se existem arquivos .in
if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo "Nenhum arquivo .in encontrado no diretório atual."
    exit 1
fi

# Mostra a lista de arquivos que serão processados
echo "-------------------------------------"
echo "Arquivos .in encontrados no diretório:"
echo "-------------------------------------"
for input_file in "${INPUT_FILES[@]}"; do
    echo "- $(basename "$input_file")"
done
echo "-------------------------------------"
echo "Total de arquivos a processar: ${#INPUT_FILES[@]}"
echo "-------------------------------------"
