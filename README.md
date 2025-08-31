# Base de dados Óxido de Zinco para NNPs

Este repositório contém a geração de uma base de dados de estrutura eletrônica e forças atômicas do óxido de zinco (ZnO) utilizando o Quantum ESPRESSO.

O objetivo será de construir uma base de dados de estruturas de ZnO com diferentes geometrias e configurações, obtidas via cálculos ab initio (DFT), com:

- **Energias totais**
- **Forças atômicas**
- **Geometrias distorcidas e/ou próximas ao equilíbrio**

Esses dados serão utilizados para **treinamento e validação de potenciais de campo de força baseados em redes neurais**, como os usados pelo [ænet](https://ann.atomistic.net/), para futuras aplicações em simulações de Dinâmica Molecular.

---

## 📁 Estrutura do Repositório
__ZnO_database/__

├── `ZnO_system_preparation/`  
   Preparação do sistema de ZnO: relaxação, testes de convergência, benchmarks  

├── `data/`  
Inputs e outputs principais das simulações DFT.  

├── `scripts/`                
 Scripts Python auxiliares.  

├── `templates_QE/`           
 Templates de input para Quantum ESPRESSO  

├── `pseudos/`  
 Pseudopotenciais utilizados (.UPF)  

├── `Notebooks/`  
Notebooks de análise e visualização dos resultados  

├── `markdown_notes/`         
 Documentações auxiliares e anotações em Markdown  

└── __README.md__  

---

## Fluxo de Trabalho

1. **Preparação do sistema (`ZnO_system_preparation/`)** - _Finalizado_
   - Relaxações iniciais
   - Distúrbios aleatórios
   - Geração de geometrias iniciais

2. **Geração de supercélulas**
   - ±10% _strain anisotrópico_ da configuração de equilíbrio para supercélulas de 1x1x1 a 3x3x3 em `anistropic_strain-Aug22`.
   
   

3. **Cálculos DFT**
   - Inputs e Outputs organizados em `data/`
   - Execução via script em `scripts/run-qe.py`

4. **Pós-processamento**
   - Extração de forças, energias e posições finais
   - Conversão para formatos aceitos por ænet

---

##  Dependências e Ferramentas

- [Quantum ESPRESSO](https://www.quantum-espresso.org/)
- [ASE (Atomic Simulation Environment)](https://wiki.fysik.dtu.dk/ase/)
- [ænet](https://ann.atomistic.net/)
- Python ≥ 3.8 com bibliotecas:
  - `numpy`, `ase`, `matplotlib`, `pandas`, etc.

---

## Autor

João Victor Campos Costa
Mestrando em Física — Universidade Federal de Minas Gerais
Contato: `victorjvc2020@ufmg.br` | `victorjvc2020@gmail.com`



