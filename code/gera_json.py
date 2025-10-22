import random
from typing import Dict, Any, Optional
import os
from pathlib import Path
import json

PATH_DATA = Path('../data')
INPUT_SEQUENCES = PATH_DATA / 'full_dataset_plasmodium'

def select_random_fasta_files() -> list[str]:
    """
    Seleciona aleatoriamente um número entre 1 e 5 de arquivos FASTA de um conjunto predefinido.
    Retorna uma lista de nomes de arquivos.
    """
    fasta_files = [file for file in os.listdir(INPUT_SEQUENCES) if file.endswith('.fasta')]
    num_files = random.randint(2, len(fasta_files))
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

def gera_parametros_aleatorios():
    """
    Gera um arquivo JSON com parâmetros aleatórios para alinhamento de sequências.
    O arquivo é salvo como 'config.json'.
    """

    params = {"algoritmo": "clustalw"}
    params["parametros"] = random_clustalw_align_params()
    params["tree_format"] = random.choice(["nexus", "newick"])
    params["entradas"] = select_random_fasta_files()

    with open("config.json", "w", encoding="utf-8") as arquivo:
        json.dump(params, arquivo, indent=4, ensure_ascii=False)


if __name__ == "__main__":
    # gera parâmetros aleatórios com seed para reprodutibilidade
    gera_parametros_aleatorios()