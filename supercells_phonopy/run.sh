#!/bin/bash

# Nome do executável do Quantum Espresso (ajuste conforme necessário)
QE_EXEC="../../../bin/pw.x"

# Número de processadores por job
NP=8

# Loop através de todos os arquivos .in
for input_file in ZnO-*.in
do
    # Remove a extensão .in para criar o nome do arquivo de output
    base_name="${input_file%.in}"
    
    # Comando para executar o Quantum Espresso
    cmd="mpirun -np $NP $QE_EXEC -in $input_file > ${base_name}.out"
    
    # Executa o comando nohup em background
    echo "Iniciando cálculo para $input_file..."
    bash -c "$cmd" > ${base_name}.out 2>&1 &
done

echo "Todos os jobs foram iniciados em background."
echo "Verifique os arquivos .out e .nohup.out para o progresso."
