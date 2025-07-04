#!/bin/bash

# =================== CONFIGURAÇÕES =======================
export OMP_NUM_THREADS=2                    # 2 threads por processo MPI (uso de SMT)
QE_EXEC=../../../../bin/pw.x                # Caminho do executável do QE
MAX_JOBS=2                                  # Número de simulações simultâneas
NP_PER_JOB=8                                # MPI por simulação (np)
NPOOL=2                                     # npool do QE (ajuste conforme seus k-points)
# ==========================================================

run_qe() {
    input="$1"
    base="${input%.in}"

    echo "Iniciando simulação: $input"
    mpirun -np "$NP_PER_JOB" \
        --map-by core:PE="$OMP_NUM_THREADS" \
        --bind-to core \
        "$QE_EXEC" -in "$input" -npool "$NPOOL" > "${base}.out" 2> "${base}.err"
    echo "Finalizada simulação: $input"
}

job_count=0

# Loop por todos os arquivos ZnO-*.in em ordem alfanumérica
for input in ZnO-*.in; do
    run_qe "$input" &
    ((job_count++))

    if [[ $job_count -ge $MAX_JOBS ]]; then
        wait
        job_count=0
    fi
done

wait
echo "Todas as simulações foram concluídas."
