# Base de dados Óxido de Zinco para NNPs

Este repositório contém a geração de uma base de dados de estrutura eletrônica e forças atômicas do óxido de zinco (ZnO) utilizando o Quantum ESPRESSO.  

O objetivo será de construir uma base de dados de estruturas de ZnO com diferentes geometrias e configurações, obtidas via cálculos ab initio (DFT), com:

- **Energias totais**
- **Forças atômicas**
- **Geometrias distorcidas ou relaxadas**

Esses dados serão utilizados para **treinamento e validação de potenciais de campo de força baseados em redes neurais**, como os usados pelo [ænet](https://ann.atomistic.net/), para futuras aplicações em simulações de Dinâmica Molecular.

---

## 📁 Estrutura do Repositório
__ZnO_database/__  
├── `ZnO_system_preparation/` # Pré-processamento do sistema (relaxações, distorções, etc.)  
├── `data/`                   # Inputs e outputs principais das simulações DFT  
├── `SCF-RND025-COMPARE/`     # Teste inputs com deslocamentos aleatórios com $ \text{ibrav=4}$ e $ \text{ibrav=6}$  
├── `scripts/`                # Scripts Python para geração, execução e ETL.  
├── `templates_QE/`           # Templates de input para Quantum ESPRESSO  
├── `pseudos/`                # Pseudopotenciais utilizados (.UPF)  
├── `supercells_phonopy/`     # Supercélulas geradas com o Phonopy  
├── `supercells-ASE/`         # Supercélulas geradas com o ASE  
├── `post-processing/`        # Scripts e ferramentas de pós-processamento  
├── `Notebooks/`              # Notebooks de análise e visualização dos resultados  
├── `markdown_notes/`         # Documentações auxiliares e anotações em Markdown  
└── __README.md__   

---

## Fluxo de Trabalho

1. **Preparação do sistema (`ZnO_system_preparation/`)** - _Finalizado_
   - Relaxações iniciais
   - Distúrbios aleatórios
   - Geração de geometrias iniciais

2. **Geração de supercélulas**
   - Via ASE (`supercells-ASE/`)
   - Via Phonopy (`supercells_phonopy/`)

3. **Cálculos DFT**
   - Entradas organizadas em `data/`
   - Execução via scripts em `scripts/`

4. **Pós-processamento**
   - Extração de forças, energias e posições finais
   - Conversão para formatos aceitos por ænet

5. **Treinamento do potencial**
   $$\vec{x_i} = f_a^i \big(\mathbf{W}_i \vec{x}_{i-1} + \vec{b}_i \big) \\ \text{Node $i$ - ANN}$$ 
   
   - Número de layers
   - Nodes per layer
   - Função de ativação.


---

##  Dependências e Ferramentas

- [Quantum ESPRESSO](https://www.quantum-espresso.org/)
- [ASE (Atomic Simulation Environment)](https://wiki.fysik.dtu.dk/ase/)
- [Phonopy](https://phonopy.github.io/phonopy/)
- [ænet](https://ann.atomistic.net/)
- Python ≥ 3.8 com bibliotecas:
  - `numpy`, `ase`, `matplotlib`, `pandas`, etc.

---


## Autor

João Victor Campos Costa 
Mestrando em Física — Universidade Federal de Minas Gerais 
Contato: `victorjvc2020@ufmg.br` | `victorjvc2020@gmail.com`



