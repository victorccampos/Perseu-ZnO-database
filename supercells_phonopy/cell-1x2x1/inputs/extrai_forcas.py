import argparse
import glob
import os

def extrair_forcas_qe(caminho_arquivo):
    bloco_forcas = []
    capturando = False

    with open(caminho_arquivo, 'r') as f:
        for linha in f:
            if "Forces acting on atoms" in linha:
                capturando = True
                continue
            if capturando:
                if linha.strip() == "" or "contrib" in linha:
                    break
                bloco_forcas.append(linha.rstrip())

    return bloco_forcas


def main():
    parser = argparse.ArgumentParser(description='Extrai blocos de forças dos arquivos de saída do QE.')
    parser.add_argument('-in', '--input', nargs='+', required=True, help='Arquivos .out de entrada (ex: ZnO-*.out)')
    args = parser.parse_args()

    arquivos = []
    for pattern in args.input:
        arquivos.extend(glob.glob(pattern))

    if not arquivos:
        print("Nenhum arquivo encontrado com os padrões fornecidos.")
        return

    with open("forcas_extraidas.txt", "w") as f_out:
        for arq in sorted(arquivos):
            f_out.write(f"{os.path.basename(arq)}\n")
            bloco = extrair_forcas_qe(arq)
            if bloco:
                f_out.write("\n".join(bloco))
            else:
                f_out.write("[Nenhum bloco de forças encontrado]")
            f_out.write("\n\n")  # Separador entre arquivos


if __name__ == "__main__":
    main()
