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

# Pede confirmação antes de continuar
#read -p "Deseja iniciar a execução? (s/n) " -n 1 -r
#echo
#if [[ ! $REPLY =~ ^[Ss]$ ]]; then
#    echo "Execução cancelada pelo usuário."
#    exit 0
#fi

# Número de processos MPI
NP=16
# Número de pools
NPOOL=2
# Opções de mapeamento
MAP_OPTIONS="--map-by numa:pe=2"

# Executa cada cálculo sequencialmente
for input_file in "${INPUT_FILES[@]}"; do
    # Remove o caminho do diretório, ficando apenas com o nome do arquivo
    file_name=$(basename "$input_file")
    
    # Nome base sem extensão
    base_name="${file_name%.*}"
    
    # Arquivos de saída e erro
    output_file="${base_name}.out"
    error_file="${base_name}.err"
    
    echo "Iniciando cálculo para ${file_name}..."
    echo "Saída será gravada em ${output_file}"
    echo "Erros serão gravados em ${error_file}"
    
    # Comando de execução
    mpirun -np $NP $MAP_OPTIONS ../../bin/pw.x -in "$file_name" -npool $NPOOL > "$output_file" 2> "$error_file"
    
    # Aguarda a conclusão do trabalho anterior antes de iniciar o próximo
    wait
    
    echo "Cálculo para ${file_name} concluído."
    echo "-------------------------------------"
done

echo "Todos os cálculos foram concluídos."
