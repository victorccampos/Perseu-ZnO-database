import os
import argparse
import numpy as np

#==============================================================================#
#                  
#                  celldm(1) = 6.178821408099141 a
#                  celldm(3) = 1.614358356153010 c/a  
#   
# 

BASE_POSITIONS_CRYSTAL = """
Zn            0.6666666667        0.3333333333        0.5006335854
Zn            0.3333333333        0.6666666667        0.0006335854
O             0.6666666667        0.3333333333        0.8796764146
O             0.3333333333        0.6666666667        0.3796764146
"""

celldm1 = 6.178821408099141 # Bohr Radius
celldm3 = 1.614358356153010 # Adimensional

# -10% a 10% dos parâmetros de rede
percent_range = np.arange(start=-0.10, stop=0.12, step=0.02)  # [start, stop)

a_range = celldm1 * (1 + percent_range)  # celldm1 * [0.9  0.92 0.94 0.96 0.98 1.   1.02 1.04 1.06 1.08 1.1 ]
ca_ratio_range = celldm3 * (1 + percent_range) # celldm1 * [0.9  0.92 0.94 0.96 0.98 1.   1.02 1.04 1.06 1.08 1.1 ]


#==============================================================================#
NOISE_LEVEL = 0.15

def create_qe_input(a_val: float, ca_ratio: float,
                     template_str: str, noise_flag: bool = False, noise_level = NOISE_LEVEL) -> str:
    """Gera o conteúdo do arquivo de input do Quantum Espresso com ruído.
    a_val: valor do parâmetro de rede a entre ±10% * CELLDM1
    ca_ratio: valor da razão c/a entre ±10% * CELLDM3
    template_str: input do QE com placeholders
            - __CELLDM1__
            - __CELLDM3__
            - __ATOMIC_POSITIONS__ (crystal coords)
    noise_flag: Adição ou não de ruído em __ATOMIC_POSITIONS__
    noise_level: deslocamento em ± noise_level em BASE_POSITIONS (crystal coords)
    """

    noisy_positions_list = []
    
    for line in BASE_POSITIONS.strip().split('\n'):
        parts = line.split()
        symbol = parts[0]
        coords = np.array([float(c) for c in parts[1:]])
        
        # Gera ruído no intervalo [-noise_level, +noise_level]
        noise = np.random.uniform(-noise_level, noise_level, 3) if noise_flag else 0.
        
        new_coords = coords + noise
        noisy_positions_list.append(f"{symbol:2s}  {new_coords[0]:.12f}  {new_coords[1]:.12f}  {new_coords[2]:.12f}")

    noisy_positions_str = "\n".join(noisy_positions_list)

    # Substitui as flags nos placeholders do template
    content = template_str.replace('__CELLDM1__', f"{a_val:.12f}")
    content = content.replace('__CELLDM3__', f"{ca_ratio:.12f}")
    content = content.replace('__ATOMIC_POSITIONS__', noisy_positions_str)

    return content

def main():
    """Função principal para gerar a estrutura de diretórios e os arquivos."""
    
    parser = argparse.ArgumentParser(
        description="Gera diretório com inputs do QE com variações de parâmetros")
    
    # Folder Flag
    parser.add_argument("-d", "--directory",
                        default="SFCs.in",
                        help="Nome diretório onde serão salvos I/O do QE")

    # Noise Flag
    parser.add_argument("-n", "--noise",
                        action="store_true", # Se passado como argumento -> True
                        help="Adiciona ruído aleatório às pos. atômicas (Padrão: False)")


    parser.add_argument("--num_configs", type=int, default=5, help="Número de configurações aleatórias por (a, c/a)")

    # Template de input - com simetria (ibrav = 4) ou sem simetria (ibrav=0 + CELL_PARAMETERS)
    parser.add_argument("-t", "--template", 
                        default="template.in", 
                        help="Template de arquivo '.in' do QE (sym/no-sym).\
                              Default: 'template.in' (ibrav = 4)")
    # Argumentos do parser
    args = parser.parse_args()


    try:
        # template_com_simetria = 'template.in' # arquivo ibrav = 4.
        # template_sem_simetria = 'template_no-sym.in'
        
        with open(args.template, 'r') as f:
            template_content = f.read()
    
    except FileNotFoundError:
        print(f"Erro: Arquivo {args.template} não encontrado. Crie-o no mesmo diretório.")
        return

    print("Iniciando a geração de arquivos...")
    print(f"Modo com ruído: {'On' if args.noise else 'Off'}")

    
    os.makedirs(args.directory, exist_ok=True)

    for config in range(1, args.num_configs + 1):  # Gera configurações de 1 a N
        for a in a_range:                           
            for ca in ca_ratio_range:               
                
                # string placeholders for files
                a_str = f"{a:.2f}"
                ca_str = f"{ca:.2f}"


                # Gera o conteúdo do input
                input_content = create_qe_input(a, ca, template_content, noise_flag=args.noise)

                input_filename: str = f"scf-{a_str}-{ca_str}-{NOISE_LEVEL}-{config}.in" if args.noise else f"scf-{a_str}-{ca_str}-.in"
                
                # Escreve o arquivo de input
                file_path = os.path.join(args.directory,input_filename )
                
                with open(file_path, 'w') as f:
                    f.write(input_content)

                print(f" Arquivo criado: {file_path}")

        print("\nProcesso concluído com sucesso!")

if __name__ == "__main__":
    main()
