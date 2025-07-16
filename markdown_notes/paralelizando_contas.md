Com base no seu hardware (AMD Opteron 6282 SE, 4 sockets, 8 núcleos por socket, 2 threads por núcleo, total de 64 CPUs lógicas) e nos cálculos de ZnO (de 4 a 108 átomos), aqui estão as **diretivas de paralelismo otimizadas** para o Quantum ESPRESSO (QE), considerando **MPI + OpenMP híbrido**:

---

### **1. Estratégia Geral**
- **MPI**: Distribui processos entre **sockets/NUMA nodes** (evita comunicação lenta entre nós).  
- **OpenMP**: Usa threads dentro de um socket para acelerar operações locais (FFTs, álgebra linear).  
- **Prioridade**:  
  - Sistemas pequenos (≤ 32 átomos): Foco em paralelização de **k-points (`-nk`)** e **FFTs (`-ntg`)**.  
  - Sistemas grandes (> 32 átomos): Ativar paralelização de **bandas (`-nb`)** e **álgebra linear (`-nd`)**.  

---

### **2. Configurações Recomendadas por Faixa de Tamanho**
#### **(A) Pequenas supercélulas (4–24 átomos)**
- **Baixo custo computacional**: Paralelize principalmente **k-points** e **FFTs**.  
- **Comando sugerido**:  
  ```bash
  mpirun -np 16 pw.x -nk 4 -ntg 4 -i input.in
  ```
  - **`-np 16`**: 16 processos MPI (1 por núcleo físico, 2 sockets).  
  - **`-nk 4`**: 4 pools de k-points (divide 216 k-points em 4 grupos de 54).  
  - **`-ntg 4`**: 4 task groups para FFTs (cada um com 4 processos MPI).  
  - **OpenMP**: Desativado (ou 2 threads por processo se habilitado).  

#### **(B) Médias supercélulas (32–72 átomos)**
- **Cálculos mais pesados**: Adicione paralelização de **bandas** e **álgebra linear**.  
- **Comando sugerido**:  
  ```bash
  mpirun -np 32 pw.x -nk 8 -ntg 4 -nb 2 -nd 16 -i input.in
  ```
  - **`-np 32`**: 32 processos MPI (8 por socket).  
  - **`-nk 8`**: 8 pools de k-points (27 k-points por pool).  
  - **`-nb 2`**: 2 grupos de bandas.  
  - **`-nd 16`**: 16 processos para álgebra linear (ScaLAPACK/ELPA).  

#### **(C) Grandes supercélulas (108 átomos)**
- **Máximo paralelismo**: Use todos os níveis (k-points, bands, FFTs, álgebra linear).  
- **Comando sugerido**:  
  ```bash
  mpirun -np 64 pw.x -nk 8 -ntg 8 -nb 4 -nd 36 -i input.in
  ```
  - **`-np 64`**: Todos os núcleos lógicos (2 threads por núcleo físico).  
  - **`-nd 36`**: Grid 6×6 para ScaLAPACK (maior eficiência em matrizes grandes).  

---

### **3. Ajustes para NUMA e OpenMP**
- **Bind de processos a sockets**: Use `--map-by socket` (MPI) e `--bind-to socket` para evitar saltos entre NUMA nodes.  
  ```bash
  mpirun --map-by socket --bind-to socket -np 32 ...
  ```
- **OpenMP**: Se habilitado, defina:  
  ```bash
  export OMP_NUM_THREADS=2  # 2 threads por processo MPI (1 por núcleo físico)
  export OMP_PROC_BIND=true
  ```

---

### **4. Exemplo Completo (108 átomos)**
```bash
# Configuração híbrida MPI+OpenMP
export OMP_NUM_THREADS=2
mpirun --map-by socket --bind-to socket -np 64 pw.x -nk 8 -ntg 8 -nb 4 -nd 36 -i ZnO_108atoms.in
```
- **Justificativa**:  
  - **`-nk 8`**: Divide 216 k-points em 8 pools (27 k-points cada).  
  - **`-ntg 8`**: 8 task groups para FFTs (8 processos MPI cada).  
  - **`-nb 4`**: 4 grupos de bandas (útil para híbridos PBE0).  
  - **`-nd 36`**: Grid 6×6 para ScaLAPACK (evita gargalo na diagonalização).  

---

### **5. Dicas Adicionais**
1. **Teste pequeno primeiro**:  
   ```bash
   mpirun -np 4 pw.x -nk 2 -ntg 2 -i ZnO_4atoms.in
   ```
   Verifique se o cálculo converge e ajuste `mixing_beta` se necessário.  

2. **Monitoramento**: Use `top` ou `htop` para verificar carga CPU e balanceamento.  

3. **I/O**: Para sistemas grandes, use `-ndiag 1` se E/S for um gargalo.  

4. **Compilação**: Certifique-se de que o QE foi compilado com:  
   ```bash
   --enable-openmp --with-scalapack
   ```

---

### **Resumo das Flags por Tamanho**
| **Átomos** | **`-np`** | **`-nk`** | **`-ntg`** | **`-nb`** | **`-nd`** | **OpenMP** |  
|------------|----------|----------|-----------|----------|----------|------------|  
| 4–24       | 16       | 4        | 4         | 1        | 1        | Não        |  
| 32–72      | 32       | 8        | 4         | 2        | 16       | Opcional   |  
| 108        | 64       | 8        | 8         | 4        | 36       | Sim (2T)   |  

Se precisar de ajustes finos para um caso específico, me avise!