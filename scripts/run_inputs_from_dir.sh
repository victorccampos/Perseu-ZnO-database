#!/bin/bash

# Configurações padrão
INPUT_DIR="."
OUTPUT_DIR="."
QE_EXEC="pw.x"
NP=16
MPI_OPTIONS="--map-by numa:pe=2"  # Otimização para sistemas NUMA
QE_OPTIONS="-npool 4"             # Paralelização interna do QE
LOG_FILE="qe_run.log"

# Processa argumentos
while [[ $# -gt 0 ]]; do
    case "$1" in
        -in)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Erro: Opção inválida $1"
            exit 1
            ;;
    esac
done

# Verifica diretórios
mkdir -p "$OUTPUT_DIR"
echo "--- Registro de Execução $(date) ---" > "$OUTPUT_DIR/$LOG_FILE"
echo "Input: $INPUT_DIR | Output: $OUTPUT_DIR" >> "$OUTPUT_DIR/$LOG_FILE"

# Contadores
SUCCESS=0
FAILED=0

for input_file in "$INPUT_DIR"/ZnO-*.in; do
    base_name=$(basename "${input_file%.in}")
    output_file="$OUTPUT_DIR/${base_name}.out"

    echo "Executando: $base_name..."
    echo "$(date) - Início: $base_name" >> "$OUTPUT_DIR/$LOG_FILE"

    # COMANDO OTIMIZADO (nohup + MPI + QE options)
    if nohup mpirun -np $NP $MPI_OPTIONS $QE_EXEC $QE_OPTIONS -in "$input_file" > "$output_file" 2>&1; then
        status="SUCESSO"
        ((SUCCESS++))
    else
        status="FALHA"
        ((FAILED++))
    fi

    echo "$(date) - $status: $base_name" >> "$OUTPUT_DIR/$LOG_FILE"
done

# Resumo final
echo "--------------------------------" >> "$OUTPUT_DIR/$LOG_FILE"
echo "Total: $((SUCCESS + FAILED)) jobs" >> "$OUTPUT_DIR/$LOG_FILE"
echo "Sucessos: $SUCCESS | Falhas: $FAILED" >> "$OUTPUT_DIR/$LOG_FILE"
echo "Log salvo em: $OUTPUT_DIR/$LOG_FILE"

# Resultado no terminal
echo -e "\nConcluído! Resultados:"
echo "• Arquivos de saída em: $OUTPUT_DIR (ZnO-*.out)"
echo "• Log detalhado em: $OUTPUT_DIR/$LOG_FILE"
echo "• Jobs executados: $SUCCESS sucesso(s), $FAILED falha(s)"
