import random
from typing import Dict, Any, Optional
import os
from pathlib import Path
import json
import numpy as np

PATH_DATA = Path('../data')
INPUT_SEQUENCES = PATH_DATA / 'full_dataset_plasmodium'

def sorteio_exponencial_inteiro(limite=400, scale=80):
    """
    Sorteia um número inteiro entre 0 e `limite`,
    com maior chance de valores próximos de zero,
    usando uma distribuição exponencial.

    Parâmetros:
        limite (int): valor máximo possível (default = 400)
        scale (float): parâmetro de escala da distribuição exponencial 
                       (valores menores concentram mais em 0)

    Retorna:
        int: número sorteado
    """
    valor = np.random.exponential(scale=scale)
    valor = int(round(valor))
    return min(valor, limite)

def select_random_fasta_files() -> list[str]:
    """
    Seleciona aleatoriamente um número entre 1 e 5 de arquivos FASTA de um conjunto predefinido.
    Retorna uma lista de nomes de arquivos.
    """
    fasta_files = [file for file in os.listdir(INPUT_SEQUENCES) if file.endswith('.fasta')]
    num_files = sorteio_exponencial_inteiro(limite=50, scale=80)
    return random.sample(fasta_files, num_files)

def random_clustalw_align_params(seed: Optional[int] = None) -> Dict[str, Any]:
    """
    Gera um dicionário com parâmetros aleatórios apropriados para um `clustalw2 -ALIGN`.
    Retorna algo no formato: {'INFILE': '...', 'TYPE': 'PROTEIN', 'ALIGN': True, ...}
    """
    if seed is not None:
        random.seed(seed)

    # Tipo das sequências
    seq_type = random.choice(["PROTEIN"]) # "DNA"

    params = {
        "-TYPE": seq_type,
        "-ALIGN": True,                      # verb ALIGN
        "-OUTPUT": random.choice(["CLUSTAL"]), # "FASTA", "PHYLIP", "NEXUS"
        "-OUTORDER": random.choice(["ALIGNED", "INPUT"]),
        "-GAPOPEN": round(random.uniform(5.0, 20.0), 1),   # ex.: 10.0
        "-GAPEXT": round(random.uniform(0.0, 1.0), 2),     # ex.: 0.20
        # "-NUMITER": random.randint(0, 8),
        # "-MAXSEQLEN": random.randint(1000, 20000),
        "-QUIET": random.choice([True, False]),
    }

    return params

def random_probcons_align_params(seed: Optional[int] = None) -> Dict[str, Any]:
    """
    Gera um dicionário com parâmetros aleatórios apropriados para um `clustalw2 -ALIGN`.
    Retorna algo no formato: {'INFILE': '...', 'TYPE': 'PROTEIN', 'ALIGN': True, ...}
    """
    if seed is not None:
        random.seed(seed)

    params = {
        "-clustalw": True,
        "-c": random.randint(0, 5),
        "-ir": random.randint(0, 1000),
        "-pre": random.randint(0, 20)
    }

    return params

def gera_parametros_aleatorios(algoritmo: str = "clustalw"):
    """
    Gera um arquivo JSON com parâmetros aleatórios para alinhamento de sequências.
    O arquivo é salvo como 'config.json'.
    """

    if algoritmo.lower() == 'clustalw':
        params = {"algoritmo": "clustalw"}
        params["parametros"] = random_clustalw_align_params()
    elif algoritmo.lower() == 'probcons':
        params = {"algoritmo": "probcons"}
        params["parametros"] = random_probcons_align_params()


    params["tree_format"] = random.choice(["nexus", "newick"])
    params["entradas"] = select_random_fasta_files()

    with open("config.json", "w", encoding="utf-8") as arquivo:
        json.dump(params, arquivo, indent=4, ensure_ascii=False)


if __name__ == "__main__":
    # gera parâmetros aleatórios com seed para reprodutibilidade
    gera_parametros_aleatorios('probcons')
