# Gerar Inputs
python3 gerar_inputs.py --template template.in --directory scfs_rnd_015.in --noise


# Rodar inputs
1. nohup python3 run_all_inputs.py -d SCFs_rnd.in/ > log_rnd.txt 2>&1 &
2. disown


# Executei assim na Barbara e deu certo 
 1122  nohup python3 run_parallel_inputs.py -d scfs_rnd025.in/ > log_rnd025.out 2>&1 &
 1123  disown

# Identificar o PID do processo
(base) jvc@perseu:~/QEspresso7.2/ZnO_database_python/scfs_rnd025-nosym.out$ ps aux | grep "run_parallel_inputs.py"
jvc       822246  0.0  0.0 612808 14336 ?        Sl   10:00   0:00 python3 run_parallel_inputs.py -d scfs_rnd025-nosym.in/
jvc       823124  0.0  0.0   9312  2048 pts/66   S+   10:34   0:00 grep --color=auto run_parallel_inputs.py

Matar o processo com kill 822246 e depois killall pw.x.

# Processo individual
mpirun -np 8 ../../bin/pw.x -in scf-rnd-ibrav0.in > scf-rnd-ibrav0.out 2> scf-rnd-ibrav0.err &

# awk para pw.x
awk '/kinetic-energy/{ecut=$4}
     /!.*total/{etot=$5}
     /Total force/{totfor=$4}
     /total.*stress/{print ecut, etot, totfor, $6}' *out
