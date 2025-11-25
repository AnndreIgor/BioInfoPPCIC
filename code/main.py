# %%
from pathlib import Path
import time
import numpy as np
import os
from itertools import combinations
import sqlite3
import json
import datetime as dt
import random
import shutil

# %%
import parsl
from parsl import python_app, bash_app
from parsl.config import Config
from parsl.executors import ThreadPoolExecutor

# %%
from Bio import Phylo

# %%
from gera_json import gera_parametros_aleatorios
from analise_fasta import analisar_fasta
from system_info import system_summary

# %%
gera_parametros_aleatorios(random.choice(["clustalw", "probcons"]))

# %%
N = int(os.getenv("NSLOTS", os.cpu_count()))

# %%
D1 = Path('../data')
PROVENANCE = D1 / 'provenance'
INPUT_SEQUENCES = D1 / 'full_dataset_plasmodium'

DATA = Path(f'../data_{dt.datetime.now().strftime("%f")}')
SUBARVORES_POSSIVEIS = DATA / 'subtrees'
ARVORES_FILOGENETICAS = DATA / 'tree'
TMP = DATA / 'tmp'
DB_SIMILARIDADES = DATA / Path("similaridades.db")

DATA.mkdir(exist_ok=True)
PROVENANCE.mkdir(exist_ok=True)
SUBARVORES_POSSIVEIS.mkdir(exist_ok=True)
ARVORES_FILOGENETICAS.mkdir(exist_ok=True)
TMP.mkdir(exist_ok=True)

# %%
parsl.clear()
parsl.load(
    Config(
        executors=[ThreadPoolExecutor(max_threads=N)]
    )
)

# %%
def limpa_diretorio(diretorio: Path):
    for arquivo in diretorio.iterdir():
        if arquivo.is_file():
            arquivo.unlink()

for diretorio in [SUBARVORES_POSSIVEIS, ARVORES_FILOGENETICAS, TMP]:
    limpa_diretorio(diretorio)


DB_SIMILARIDADES.unlink(missing_ok=True)

# %%
def verifica_caracteres_validos(file: Path, valid_characters: str = 'ACDEFGHIKLMNPQRSTVWY'):
    # Validar se as sequências contêm apenas caracteres válidos
    valid_characters = set(valid_characters)

    try:
        with open(file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    continue

                sequence = line.strip().upper()
                if not set(sequence).issubset(valid_characters):
                    print(sequence)
                    return False
    except FileNotFoundError:
        print(f"Arquivo {file} não encontrado.")
        return False
    
    return True

def duplicate_names(file: Path):
    from Bio import SeqIO

    name_count = {}
    records = []
    duplicate = False
    
    try:
        for record in SeqIO.parse(file, 'fasta'):
            name_count[record.id] = name_count.get(record.id, 0) + 1

            if name_count[record.id] == 1:
                records.append(record)
            else:
                duplicate = True
        
        if duplicate:
            SeqIO.write(records, file, "fasta")
                
    except FileNotFoundError:
        print(f"O arquivo '{file}' não foi encontrado.")

@python_app
def valida_sequencias(files: list[Path]):
    for file in files:
        file = INPUT_SEQUENCES / file
        if not verifica_caracteres_validos(file):
            file.unlink()
            print(f"Arquivo '{file}' removido devido a caracteres inválidos.")
        else:
            duplicate_names(file)

# %%
def analisar_sequencias_fasta(config_file: Path) -> list:

    with open(config_file, "r", encoding="utf-8") as arquivo:
        par = json.load(arquivo)

    resultado = []

    for file in par['entradas']:
        sequencia_fasta = INPUT_SEQUENCES / file
        resultado.append(analisar_fasta(sequencia_fasta))

    par['entradas'] = resultado

    with open(DATA / "temp.json", "w", encoding="utf-8") as arquivo:
        json.dump(par, arquivo, indent=4, ensure_ascii=False)

# %%
# Aqui tem que entrar o json com os parametros e algoritmo de alinhamento
def gera_comando_alinhamento(algoritmo: str, parametros: dict, infile):
    infile = INPUT_SEQUENCES / infile

    if algoritmo.lower() == 'clustalw':
        command_parameters = ""
        for key, value in parametros.items():
            if isinstance(value, bool) and value:
                command_parameters += key
            else:
                command_parameters += f"{key}={value}"
            command_parameters += ' '
        command_parameters = command_parameters.strip()

        command_files = [f"-INFILE={infile}"]
        command_files.append(f"-OUTFILE={TMP / infile.with_suffix('.aln').name}")
        
        return  algoritmo + " " + command_parameters + " " + " ".join(command_files)
    elif algoritmo.lower() == 'probcons':
        command_parameters = ""
        for key, value in parametros.items():
            if isinstance(value, bool) and value:
                command_parameters += key
            else:
                command_parameters += f"{key} {value}"
            command_parameters += ' '
        command_parameters = command_parameters.strip()

        command_files = f"{infile} >  {TMP / infile.with_suffix('.aln').name}"

        return  algoritmo + " " + command_parameters + " " + command_files

@bash_app
def alinhamento_de_sequencias(algoritmo, paramentros, infile,  stdout='echo-hello.stdout', stderr='echo-hello.stderr') -> str:
    return gera_comando_alinhamento(algoritmo, paramentros, infile)

# %%
@python_app
def geracao_de_arvores_filogeneticas(file: Path, path_out_tree: Path, modelo_evolutivo: str = 'nj', output_format: str = 'nexus'):
    from Bio import AlignIO, Phylo
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

    with open(file, "r") as handle:
        alignment = AlignIO.read(handle, "clustal")
    
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment) 
    constructor = DistanceTreeConstructor()

    if modelo_evolutivo == 'nj':
        tree = constructor.nj(distance_matrix)
    elif modelo_evolutivo == 'upgma':
        tree = constructor.upgma(distance_matrix)

    Phylo.write(tree, path_out_tree / file.with_suffix(f".{output_format}").name, output_format)

# %%
@python_app
def geracao_de_subarvores_possiveis(file: Path, data_output_path: Path, data_format: str = "nexus"):
    from Bio import Phylo

    # Salva a árvore
    tree = Phylo.read(file, data_format)

    # Lista caminhos das subárvores (que posteriormente serão utilizadas para compor a matriz de subárvores)
    row_subtree = []

    for clade in tree.find_clades():
        subtree = Phylo.BaseTree.Tree(clade)
        if subtree.count_terminals() > 1:
            filepath_out = data_output_path / f'{file.stem}_{clade.name}.{data_format}'
            Phylo.write(subtree, filepath_out, data_format)        
            row_subtree.append(filepath_out)

# %%
def mapeamento_de_subarvores(path: Path, valor_preenchimento = None):
    grupos = {}
    for item in path.iterdir():
        sequencia = item.stem.split('_')[0]
        grupos.setdefault(sequencia, []).append(item)

    # Converte para lista de listas
    matriz = list(grupos.values())

    max_len = max(len(sublista) for sublista in matriz)
    
    for sublista in matriz:
        sublista += [valor_preenchimento] * (max_len - len(sublista))

    return np.array(matriz)

# %%
def get_subtree_names(file: Path, format: str = 'nexus'):
    subtree = Phylo.read(file, format)
    return {i.name for i in subtree.get_terminals()}

def geracao_dicionario_saida(resultado):
    similaridade = []
    contador = 0

    with sqlite3.connect(DB_SIMILARIDADES) as conn:
        cursor = conn.cursor()

        cursor.execute("""CREATE TABLE IF NOT EXISTS similaridades (
            SEQUENCIA_1 VARCHAR(50), 
            SEQUENCIA_2 VARCHAR(50), 
            GRAU_MAF FLOAT);""")

        cursor.execute("PRAGMA journal_mode = WAL;")
        cursor.execute("PRAGMA synchronous = NORMAL;")

        for comp in resultado:
            x, y = comp

            # desempacota diretamente chave e valor
            (chave_0, set_0), = x.items()
            (chave_1, set_1), = y.items()

            if chave_0.split('_')[0] == chave_1.split('_')[0]:
                continue

            g_maf = round(len(set_0 & set_1) / max(len(set_0), len(set_1)), 4)
            if g_maf > 0:
                similaridade.append((chave_0,chave_1,g_maf))

            contador += 1
            if contador == 10_000:
                cursor.executemany("""
                    INSERT INTO similaridades (SEQUENCIA_1, SEQUENCIA_2, GRAU_MAF)
                    VALUES (?, ?, ?)
                """, similaridade)

                contador = 0
                similaridade = []
        
        # Se tem um lote menor que 10_000
        if similaridade:
            cursor.executemany("""
                INSERT INTO similaridades (SEQUENCIA_1, SEQUENCIA_2, GRAU_MAF)
                VALUES (?, ?, ?)
            """, similaridade)

# %%
inicio_geral = time.perf_counter()

# %%
with open("config.json", "r", encoding="utf-8") as arquivo:
    par = json.load(arquivo)

analisar_sequencias_fasta("config.json")

# %%
inicio = time.perf_counter()

# Valida as sequencias .fasta
future = valida_sequencias(par['entradas'])

future.result()  # Aguarda a conclusão da validação

fim = time.perf_counter()
print(f'Executado em {fim - inicio:0.2f} segundos')

# %%
inicio = time.perf_counter()

# Alinha as sequencias .fasta
futures = []
for file in par['entradas']:
    # print(gera_comando_alinhamento(par['algoritmo'], par['parametros'], file))
    futures.append(alinhamento_de_sequencias(par['algoritmo'], par['parametros'], file))
    
for f in futures:
    f.result()  # Aguarda a conclusão de cada tarefa

fim = time.perf_counter()
print(f'Executado em {fim - inicio:0.2f} segundos')

# %%
inicio = time.perf_counter()

# Gera as árvores filogenéticas
futures = []
for file in TMP.glob('*.aln'):
    futures.append(geracao_de_arvores_filogeneticas(file, ARVORES_FILOGENETICAS, output_format=par['tree_format']))

for f in futures:
    f.result()  # Aguarda a conclusão de cada tarefa

fim = time.perf_counter()
print(f'Executado em {fim - inicio:0.2f} segundos')

# %%
inicio = time.perf_counter()

futures = []
for file in ARVORES_FILOGENETICAS.glob(f"*.{par['tree_format']}"):
    futures.append(geracao_de_subarvores_possiveis(file, SUBARVORES_POSSIVEIS, par['tree_format']))

for f in futures:
    f.result()  # Aguarda a conclusão de cada tarefa

fim = time.perf_counter()
print(f'Executado em {fim - inicio:0.2f} segundos')

# %%
inicio = time.perf_counter()

matriz_subtree  = mapeamento_de_subarvores(SUBARVORES_POSSIVEIS)

fim = time.perf_counter()
print(f'Executado em {fim - inicio:0.2f} segundos')

# %%
lista_subarvores = []
for subtree in SUBARVORES_POSSIVEIS.iterdir():
    lista_subarvores.append({subtree.name: get_subtree_names(file=subtree, format=par['tree_format'])})

resultado_listas = list(combinations(lista_subarvores, 2))

inicio = time.perf_counter()

geracao_dicionario_saida(resultado_listas)

fim = time.perf_counter()
print(f'Executado em {fim - inicio:0.2f} segundos')

# %%
fim_geral = time.perf_counter()
print(f'Executado em {fim_geral - inicio_geral:0.2f} segundos')

# %%
with open(DATA / "temp.json", "r", encoding="utf-8") as arquivo:
    par = json.load(arquivo)

par['resultado'] = {"Inicio": inicio_geral,
                    "Fim": fim_geral,
                    "num_procs": N}
        
par['host'] = system_summary()
        
with open(PROVENANCE / f"dados_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.json", "w", encoding="utf-8") as arquivo:
    json.dump(par, arquivo, indent=4, ensure_ascii=False)
        
os.remove(DATA / "temp.json")
shutil.rmtree(DATA)


