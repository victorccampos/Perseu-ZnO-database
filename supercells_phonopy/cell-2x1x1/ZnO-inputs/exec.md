nohup mpirun -np 8 ../../../../bin/pw.x -in ZnO-001.in > ZnO-001.out 2> ZnO-001.err &

__Execução Otimizada ?__
2000  export OMP_NUM_THREADS=2
2009  nohup mpirun -np 16 --map-by numa:pe=2 ../../../../bin/pw.x -in ZnO-001.in -npool 2 > ZnO-001.out 2> ZnO-001.err &

__Energias__
awk '/^!/ {print FILENAME, $5}' *.out > dados_energia.dat
__Forcas__
awk '/Total force/ {print FILENAME ": " $4 }' *.out >> dados_energia.dat 
