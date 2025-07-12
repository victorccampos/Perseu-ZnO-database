Aqui está o arquivo editado com as tabelas centralizadas:

### Tabela de Número de Átomos em Supercélulas de ZnO (Wurtzita)

A estrutura wurtzita do ZnO possui 4 átomos por célula unitária (2 Zn e 2 O). Abaixo está a tabela com o número total de átomos para todas as combinações de supercélulas de $1 \times 1 \times 1$ até $3 \times 3 \times 3$, variando $n_x$, $n_y$ e $n_z$:

<center>

| $$n_x$$ | $$n_y$$ | $$n_z$$ | Nº de Átomos |
|-------|-------|-------|--------------|
| 1     | 1     | 1     | 4            |
| 1     | 1     | 2     | 8            |
| 1     | 1     | 3     | 12           |
| 1     | 2     | 1     | 8            |
| 1     | 2     | 2     | 16           |
| 1     | 2     | 3     | 24           |
| 1     | 3     | 1     | 12           |
| 1     | 3     | 2     | 24           |
| 1     | 3     | 3     | 36           |
| 2     | 1     | 1     | 8            |
| 2     | 1     | 2     | 16           |
| 2     | 1     | 3     | 24           |
| 2     | 2     | 1     | 16           |
| 2     | 2     | 2     | 32           |
| 2     | 2     | 3     | 48           |
| 2     | 3     | 1     | 24           |
| 2     | 3     | 2     | 48           |
| 2     | 3     | 3     | 72           |
| 3     | 1     | 1     | 12           |
| 3     | 1     | 2     | 24           |
| 3     | 1     | 3     | 36           |
| 3     | 2     | 1     | 24           |
| 3     | 2     | 2     | 48           |
| 3     | 2     | 3     | 72           |
| 3     | 3     | 1     | 36           |
| 3     | 3     | 2     | 72           |
| 3     | 3     | 3     | 108          |

</center>

Cada valor corresponde ao total de átomos (Zn + O) na supercélula para a respectiva combinação de $n_x$, $n_y$ e $n_z$.

## Configurações Adaptadas ao Seu Hardware AMD Opteron

Baseando-se nas especificações do seu servidor[1], adaptei a tabela de configurações para otimizar o uso dos recursos disponíveis:

### Análise do Hardware

Seu servidor possui características importantes para simulações DFT:
- **64 CPUs totais** (32 núcleos físicos com 2 threads cada)
- **8 nós NUMA** com 8 CPUs cada
- **4 soquetes** com 8 núcleos por soquete
- **Arquitetura AMD Opteron 6282 SE**

### Tabela Adaptada para Supercélulas de ZnO

<center>

| Supercélula | Nº Átomos | Configuração Otimizada |
|-------------|-----------|------------------------|
| 1×1×1 | 4 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 1×1×2 | 8 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 1×2×1 | 8 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 2×1×1 | 8 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 1×1×3 | 12 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 1×3×1 | 12 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 3×1×1 | 12 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 1×2×2 | 16 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 2×1×2 | 16 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 2×2×1 | 16 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 1×2×3 | 24 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 1×3×2 | 24 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 2×1×3 | 24 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 2×3×1 | 24 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 3×1×2 | 24 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 3×2×1 | 24 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 2×2×2 | 32 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 1×3×3 | 36 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 3×1×3 | 36 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 3×3×1 | 36 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 2×2×3 | 48 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 2×3×2 | 48 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 3×2×2 | 48 | `mpirun -np 16 --map-by numa:pe=2 pw.x -npool 2` |
| 2×3×3 | 72 | `mpirun -np 24 --map-by numa:pe=2 pw.x -npool 3` |
| 3×2×3 | 72 | `mpirun -np 24 --map-by numa:pe=2 pw.x -npool 3` |
| 3×3×2 | 72 | `mpirun -np 24 --map-by numa:pe=2 pw.x -npool 3` |
| 3×3×3 | 108 | `mpirun -np 32 --map-by numa:pe=2 pw.x -npool 4` |

</center>

### Otimizações Específicas para Seu Hardware

#### **Limitação de Processos MPI**
- **Máximo de 32 processos MPI** para evitar sobrecarga do sistema
- Cada processo utiliza 2 threads (`pe=2`), totalizando até 64 threads
- Configuração conservadora que mantém estabilidade e performance

#### **Aproveitamento da Arquitetura NUMA**
- **8 nós NUMA** permitem distribuição eficiente de processos
- Mapeamento `--map-by numa:pe=2` otimiza localidade de memória
- Reduz latência de acesso entre processos e dados

#### **Estratégia de Pools**
- **Pools conservadores** (2-4) para manter eficiência
- Para supercélulas menores: `-npool 2`
- Para supercélulas maiores: `-npool 3-4`

### Recomendações Adicionais

#### **Variáveis de Ambiente**
Mantenha sempre:
```bash
export OMP_NUM_THREADS=2
```

#### **Monitoramento de Performance**
- Teste diferentes configurações para supercélulas críticas
- Monitore uso de CPU e memória durante execução
- Ajuste `-npool` conforme densidade de k-points

#### **Execução Simultânea**
Com 64 CPUs disponíveis, você pode executar **duas simulações simultâneas** de 16 processos cada, ou uma única simulação com até 32 processos para supercélulas maiores.

Esta configuração otimizada aproveita as características específicas do seu hardware AMD Opteron, garantindo eficiência computacional para suas simulações de ZnO dopado.

[1] https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/57775599/313f03f9-c106-47e4-913c-7e0a919d4a59/lscpu_perseu.txt