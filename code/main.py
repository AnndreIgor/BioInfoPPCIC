# %%
# Benchmark dê comparação de strings e benchmark de leitura e escrita

# %%
import logging
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import json
import time
import pandas as pd
import random
import sqlite3
import datetime as dt
import sys

from pathlib import Path

from concurrent.futures import ProcessPoolExecutor
from typing import Iterable, List

# %%
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import SeqIO, Phylo, AlignIO

# %%
from analise_fasta import analisar_fasta
from system_info import system_summary

# %%
PATH_DATA = Path('../data')
INPUT_SEQUENCES = PATH_DATA / 'full_dataset_plasmodium'
PATH_OUT = PATH_DATA / 'out'
SEQUENCIAS_ALINHADAS = PATH_OUT / 'tmp'
ARVORES_FILOGENETICAS = PATH_OUT / 'Trees'
SUBARVORES = PATH_OUT / 'Subtrees'
PROVENANCE = PATH_DATA / 'provenance'
SIMILARIDADES = PATH_OUT / 'Similaridades'
NOME_ARQUIVO = PROVENANCE

# %% [markdown]
# ### Gera os paramtros de sequencias e entradas aleatórias

# %%
from gera_json import gera_parametros_aleatorios
gera_parametros_aleatorios()

# %% [markdown]
# ### Pega os dados das sequencias

# %%
# Intensivo em leitura
def contar_sequencias(caminho_arquivo: Path) -> int:
    """
    Conta a quantidade de sequências em um arquivo fasta.
    A contagem é feita pelo número de linhas de cabeçalho.
    """
    with open(caminho_arquivo, 'r', encoding='utf-8') as f:
        conteudo = f.read()

    return conteudo.count(">")


# %% [markdown]
# ### Funções de limpeza

# %%
def clean_Trees():
    for name_file_trees in os.listdir(ARVORES_FILOGENETICAS):
        if name_file_trees != "file.gitkeep":
            os.remove(ARVORES_FILOGENETICAS / name_file_trees)

def clean_tmp():
    for name_file_tmp in os.listdir(SEQUENCIAS_ALINHADAS):
        if name_file_tmp != "file.gitkeep":
            os.remove(SEQUENCIAS_ALINHADAS / name_file_tmp)


def clean_NoPipe():
    for file_name in os.listdir(INPUT_SEQUENCES):
        if 'NoPipe' in file_name or file_name.endswith('.dnd'):
            os.remove(INPUT_SEQUENCES / file_name)


def clean_subtrees():
    for name_file in os.listdir(SUBARVORES):
        if name_file != "file.gitkeep":
            os.remove(SUBARVORES / name_file)

# %% [markdown]
# ### Funções de validação das sequências

# %%
# Verifica se todos os caracteres da sequência são válidos
def validate_fasta_protein(file_path: Path) -> bool:
    valid_chars = set('ACDEFGHIKLMNPQRSTVWY')
    if not file_path.exists():
        logging.error(f"Arquivo não encontrado: {file_path}")
        return False
    for line in file_path.read_text().splitlines():
        if line.startswith('>'):
            continue
        if not set(line.strip()).issubset(valid_chars):
            logging.warning(f"Sequência inválida em: {file_path}")
            return False
    return True

# Verifica se existe algum id duplicado na sequencia
def has_duplicate_ids(file_path: Path) -> bool:
    if not file_path.exists():
        logging.error(f"Arquivo não encontrado: {file_path}")
        return False
    seen = set()
    for record in SeqIO.parse(file_path, 'fasta'):
        if record.id in seen:
            logging.warning(f"ID duplicado encontrado: {record.id}")
            return True
        seen.add(record.id)
    return False

# Remove sequencias duplicadas
def remove_pipe(name: str, path_in_fasta: Path) -> Path:
    sequences = list(SeqIO.parse(path_in_fasta, "fasta"))
    unique_sequences = {str(seq.seq): seq for seq in sequences}
    output_file_tmp = path_in_fasta / f"{name}_NoPipe"
    SeqIO.write(unique_sequences.values(), output_file_tmp, "fasta")
    logging.info(f"Arquivo sem duplicatas salvo em: {output_file_tmp}")
    return output_file_tmp

# %% [markdown]
# ### Analisa as sequencias para posterior inserção em modelo

# %%
def analisar_sequencias_fasta(config_file: Path) -> list:

    with open(config_file, "r", encoding="utf-8") as arquivo:
        par = json.load(arquivo)

    entradas_original = par['entradas']

    resultado = []

    for file in par['entradas']:
        sequencia_fasta = INPUT_SEQUENCES / file
        resultado.append(analisar_fasta(sequencia_fasta))

    par['entradas'] = resultado

    with open(PROVENANCE / "temp.json", "w", encoding="utf-8") as arquivo:
        json.dump(par, arquivo, indent=4, ensure_ascii=False)

    return entradas_original

# %% [markdown]
# ### Função de alinhamento sequências

# %%
def alinhar_sequencia(infile: Path, outfile: Path):
    with open("config.json", "r", encoding="utf-8") as arquivo:
        par = json.load(arquivo)

    command = []

    command.append(par['algoritmo'])

    match par['algoritmo']:
        case "mafft":
            for key, value in par['parametros'].items():
                if isinstance(value, bool) and value:
                    command.append(key)
                else:
                    command.append(f"{key} {value}")

            command.append(f"{infile}")
            
            print("Executando comando:", ' '.join(command))

            with open(outfile, "w") as out:
                result = subprocess.run(
                    command,
                    stdout=out,
                    stderr=subprocess.PIPE,  # captura stderr
                    text=True                 # recebe str em vez de bytes
                )

            # exibe erro só se houver
            if result.stderr and result.stderr.strip():
                return "Erro durante execução:" + result.stderr.strip()

            
        case "clustalw":
            command.append(f"-INFILE={infile}")
            command.append(f"-OUTFILE={outfile}")

            for key, value in par['parametros'].items():
                if isinstance(value, bool) and value:
                    command.append(key)
                else:
                    command.append(f"{key}={value}")

            result = subprocess.run(command, capture_output=True, text=True)

            if result.stderr and result.stderr.strip():
                print("Erro durante execução:", result.stderr.strip(), file=sys.stderr)

        case "probcons":
            for key, value in par['parametros'].items():
                # if isinstance(par[key], list):
                #     value = random.choice(par[key])

                if isinstance(value, bool) and value:
                    command.append(key)
                else:
                    command.append(f"{key} {value}")

            command.append(f"{infile}")
                        
            with open(outfile, "w") as out:
                subprocess.run(command, stdout=out, stderr=subprocess.DEVNULL)
    
    return None

# %% [markdown]
# ### Função para calcular o número de comparações necessárias

# %%
def calcula_numero_comparacoes(m: int, n: int) -> int:
    """
    Calcula f(m, n) = C(m, 2) * n² = (m * (m - 1) // 2) * n ** 2.

    Parâmetros
    ----------
    m : int
        Inteiro não-negativo (m ≥ 0).
    n : int
        Inteiro (pode ser zero ou positivo).

    Retorno
    -------
    int
        Valor de f(m, n).
    """
    return (m * (m - 1) // 2) * n ** 2


# %% [markdown]
# ### Construção de árvores filogenéticas

# %%
# Constroi a arvore filogenética
def constructor_tree(arquivo_alinhamentos: Path, path_out_tree: Path, output_format):
    # Lê o arquivo do alinhamento de sequencias
    with arquivo_alinhamentos.open("r") as handle:
        alignment = AlignIO.read(handle, "clustal")
    

    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Escolha do modelo evolutivo nj ou upgma
    tree = DistanceTreeConstructor().nj(distance_matrix)
    Phylo.write(tree, path_out_tree, output_format) # Escreve em arquivo

# %% [markdown]
# ### Construção das subárvores

# %%
def sub_tree(path, name_subtree, data_format: str, data_output_path: Path, extension_format: str):
    # Salva a árvore
    tree = Phylo.read(path, data_format) # Lê de arquivo
    name_subtree = name_subtree.rsplit(".", 1)[0]

    # Lista caminhos das subárvores (que posteriormente serão utilizadas para compor a matriz de subárvores)
    row_subtree = []

    for clade in tree.find_clades(terminal=False):
        subtree = Phylo.BaseTree.Tree(clade)
        filepath_out = os.path.join(data_output_path,f'{name_subtree}_{clade.name}.{extension_format}')
        Phylo.write(subtree, filepath_out, data_format) # Escreve em arquivo        
        row_subtree.append(filepath_out)
            
    return row_subtree 

# %%
def preencher_matriz(matriz, valor_preenchimento, max_columns: int):
    # Preencher as células vazias com o valor de preenchimento

    for row in matriz:
        while len(row) < max_columns:
            row.append(valor_preenchimento)

    return matriz

# %% [markdown]
# ### Comparação entre subárvores

# %%
def _count_chunk(chunk: Iterable[str], set2: set) -> int:
    # função que roda em outro processo
    return sum(1 for x in chunk if x in set2)

def chunks_from_list(lst: List[str], n_chunks: int) -> List[List[str]]:
    # cria n_chunks balanceados (slice step pode também ser usado)
    k, m = divmod(len(lst), n_chunks)
    chunks = []
    start = 0
    for i in range(n_chunks):
        end = start + k + (1 if i < m else 0)
        chunks.append(lst[start:end])
        start = end
    return chunks

def grade_maf_parallel(path_1: Path, path_2: Path, data_format: str, n_workers: int = None) -> float:
    if path_1 is None or path_2 is None:
        return -1.0

    subtree_1 = Phylo.read(path_1, data_format)
    subtree_2 = Phylo.read(path_2, data_format)

    names1 = list({t.name for t in subtree_1.get_terminals() if t.name is not None})
    set2 = {t.name for t in subtree_2.get_terminals() if t.name is not None}

    if not names1 or not set2:
        return 0.0

    if n_workers is None:
        n_workers = max(1, os.cpu_count() - 1)

    # se poucos elementos, evita paralelismo
    if len(names1) < 10000 or n_workers == 1:
        common_count = sum(1 for x in names1 if x in set2)
    else:
        chunks = chunks_from_list(names1, n_workers)
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            futures = [ex.submit(_count_chunk, ch, set2) for ch in chunks]
            common_count = sum(f.result() for f in futures)

    max_possible_maf = min(len(names1), len(set2))
    return round(common_count / max_possible_maf, 2)


# %%
# Geração das subárvores possíveis
def subarvores_possiveis(dir: Path, tree_format: str) -> list[str]:
    # matriz com todas as subárvores
    matriz_subtree = []

    for name_file in os.listdir(dir):
        if (name_file != "file.gitkeep"):
            matriz_subtree.append(sub_tree(
                dir / name_file,
                name_file,
                tree_format,
                SUBARVORES,
                tree_format
            ))

    return matriz_subtree

# %%
def calcula_similaridade(max_rows: int, max_columns: int, matriz_subtree: list[list[Path]], tree_format: str, num_workers:int = 4) -> dict:
    dict_maf_database = {}

    for i in range(max_rows):
        for j in range(max_columns):
            for k in range(i + 1, max_rows):
                for l in range(max_columns):
                    g_maf = grade_maf_parallel(matriz_subtree[i][j], matriz_subtree[k][l], tree_format, num_workers)

                    if g_maf is not None and g_maf >= 0:
                        if g_maf not in dict_maf_database:
                            dict_maf_database[g_maf] = {}

                    if os.path.split(matriz_subtree[i][j])[1] not in dict_maf_database[g_maf]:
                        dict_maf_database[g_maf][os.path.split(matriz_subtree[i][j])[1]] = []

                    dict_maf_database[g_maf][os.path.split(matriz_subtree[i][j])[1]].append(os.path.split(matriz_subtree[k][l])[1])
    
    return dict_maf_database

# %% [markdown]
# ### Faz a busca das subárvores filogenéticas

# %%
# Isso tudo pode ser paralelizado
clean_NoPipe()
clean_tmp()
clean_Trees()
clean_subtrees()

lista_saida = []
if __name__ == '__main__':
    inicio = time.perf_counter()
    
    files = analisar_sequencias_fasta(config_file="config.json")

    # Contando manual para evitar contar sequências inválidas
    qtd_sequencias_validas = 0
        
    for file in files:
        # print(file)
        sequencia_fasta = INPUT_SEQUENCES / file
        sequencia_alinhada = SEQUENCIAS_ALINHADAS / os.path.split(sequencia_fasta.with_suffix(".aln"))[1]
        tree_format = "nexus"

        # Vai até a etapa de geração de árvore filogenética
        if validate_fasta_protein(sequencia_fasta):
            if has_duplicate_ids(sequencia_fasta):
                remove_pipe(sequencia_fasta)
                
            erro = alinhar_sequencia(sequencia_fasta, sequencia_alinhada)
            
            if erro is not None:
                # Se cair aqui é pq algum dos parametros está errado
                print(erro)
                sys.exit()

            constructor_tree(sequencia_alinhada, ARVORES_FILOGENETICAS / f"{file.replace('.fasta', f'.{tree_format}')}", 'nexus')

            qtd_sequencias_validas += 1

    tempo_sciphy = inicio - time.perf_counter()
    matriz_subtree = subarvores_possiveis(ARVORES_FILOGENETICAS, tree_format)
        
    max_columns = max(len(row) for row in matriz_subtree)
    max_rows = len(matriz_subtree)

    matriz_subtree = preencher_matriz(matriz_subtree, None, max_columns)

    numero_iteracoes = calcula_numero_comparacoes(max_rows, max_columns) # Cálculo do número de comparações

    num_workers = random.randint(2, os.cpu_count() - 1)
    dict_maf_database = calcula_similaridade(max_rows, max_columns, matriz_subtree, tree_format, num_workers)

    # Salva o dicionário de similaridades
    with open(SIMILARIDADES / f"similaridades_{dt.datetime.now().strftime('%Y%m%d%H%M%S')}.json", "w", encoding="utf-8") as arquivo:
        json.dump(dict_maf_database, arquivo, ensure_ascii=False, indent=4)


    with open(PROVENANCE / "temp.json", "r", encoding="utf-8") as arquivo:
        par = json.load(arquivo)

    par['resultado'] = {"qtdSequencias": qtd_sequencias_validas, 
                        "NumComparacoes": numero_iteracoes, 
                        "Inicio": inicio,
                        "Fim": time.perf_counter(),
                        "num_procs": num_workers}
        
    par['host'] = system_summary()
        
    with open(PROVENANCE / f"dados_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.json", "w", encoding="utf-8") as arquivo:
        json.dump(par, arquivo, indent=4, ensure_ascii=False)
        
    os.remove(PROVENANCE / "temp.json")

    clean_NoPipe()
    clean_tmp()
    clean_Trees()
    clean_subtrees()


