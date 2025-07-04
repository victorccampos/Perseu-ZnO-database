#!/bin/bash

export OMP_NUM_THREADS=2
QE_EXEC=../../../../bin/pw.x

# Número máximo de jobs simultâneos -- 16 núcleos (2x8)
MAX_JOBS=2

# Função que roda um input em background
run_qe() {
    input=$1
    base=${input%.in}
    echo "Iniciando $input"
    mpirun -np 8 --map-by numa:pe=2 $QE_EXEC -in "$input" -npool 2 > "${base}.out" 2> "${base}.err"
    echo "Finalizado $input"
}

# Contador de jobs em andamento
job_count=0

# Loop por todos os arquivos .in
for input in ZnO-*.in; do
    run_qe "$input" &
    ((job_count++))

    # Se já temos MAX_JOBS rodando, aguarda todos terminarem
    if [[ $job_count -ge $MAX_JOBS ]]; then
        wait
        job_count=0
    fi
done

# Espera por qualquer restante
wait
