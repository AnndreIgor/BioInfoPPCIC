# %%
# Benchmark dê comparação de strings e benchmark de leitura e escrita

# %%
import logging
import subprocess
import os
import json
import time
import datetime as dt
import sys

from pathlib import Path

# %%
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import SeqIO, Phylo, AlignIO

# %%
from analise_fasta import analisar_fasta
from system_info import system_summary

# %%
import parsl
from parsl import python_app, ThreadPoolExecutor
from concurrent.futures import wait

# %%
NUCLEOS = int(os.getenv("NSLOTS", "4"))

# %%
PATH_DATA = Path('../data')
INPUT_SEQUENCES = PATH_DATA / 'full_dataset_plasmodium'
PATH_OUT = PATH_DATA / 'out'
SEQUENCIAS_ALINHADAS = PATH_OUT / 'tmp'
ARVORES_FILOGENETICAS = PATH_OUT / 'Trees'
SUBARVORES = PATH_OUT / 'Subtrees'
PROVENANCE = PATH_DATA / 'provenance'
SIMILARIDADES = PATH_OUT / 'Similaridades'

# %%
# Criar diretórios de saída, se não existirem
os.makedirs(PATH_OUT, exist_ok=True)
os.makedirs(SEQUENCIAS_ALINHADAS, exist_ok=True)
os.makedirs(ARVORES_FILOGENETICAS, exist_ok=True)
os.makedirs(SUBARVORES, exist_ok=True)
os.makedirs(SIMILARIDADES, exist_ok=True)
os.makedirs(PROVENANCE, exist_ok=True)

# %%
# --- Configuração do Parsl ---
parsl.load(
    parsl.config.Config(
        executors=[ThreadPoolExecutor(max_threads=NUCLEOS)],
        strategy=None
    )
)

# %% [markdown]
# ### Gera os paramtros de sequencias e entradas aleatórias

# %%
from gera_json import gera_parametros_aleatorios
gera_parametros_aleatorios()

# %% [markdown]
# ### Utilitárias

# %%
def dividir_lista(lista, n_partes):
    k, m = divmod(len(lista), n_partes)
    return [lista[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n_partes)]


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
@python_app
def alinhar_sequencias(parametros: dict, files: list[Path]) -> None:
    match parametros['algoritmo']:
        case "clustalw":
            command_parameters = []

            for key, value in parametros['parametros'].items():
                if isinstance(value, bool) and value:
                    command_parameters.append(key)
                else:
                    command_parameters.append(f"{key}={value}")

            for infile in files:
                command_files = [f"-INFILE={infile}"]
                command_files.append(f"-OUTFILE={SEQUENCIAS_ALINHADAS / infile.with_suffix('.aln').name}")

                command = ["clustalw"] + command_parameters + command_files
                result = subprocess.run(command, capture_output=True, text=True)

                if result.stderr and result.stderr.strip():
                    print("Erro durante execução:", result.stderr.strip(), file=sys.stderr)

# %% [markdown]
# ### Construção de árvores filogenéticas

# %%
@python_app
def constructor_tree(arquivos_alinhamento: list[Path], output_format):
    for arquivo_alinhamentos in arquivos_alinhamento:
        with arquivo_alinhamentos.open("r") as handle:
            alignment = AlignIO.read(handle, "clustal")
            
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)

        # Escolha do modelo evolutivo nj ou upgma
        tree = DistanceTreeConstructor().nj(distance_matrix)

        # Escreve o arquivo de saida
        path_out_tree = ARVORES_FILOGENETICAS / arquivo_alinhamentos.with_suffix("." + output_format).name
        Phylo.write(tree, path_out_tree, output_format) # Escreve em arquivo

# %% [markdown]
# ### Construção das subárvores

# %%
# @python_app
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
def subarvores_possiveis(dir: Path, tree_format: str) -> list[str]:
    # matriz com todas as subárvores
    matriz_subtree = []

    for file in ARVORES_FILOGENETICAS.iterdir():
        matriz_subtree.append(
            sub_tree(
                file,
                file.name,
                tree_format,
                SUBARVORES,
                tree_format
            )
        )

    return matriz_subtree

# %%
def grade_maf(path_1: str, path_2: str, tree_format: str) -> float:
    # Verifica se os caminhos das subárvores existe
    if path_1 is None or path_2 is None:
        return None
    
    subtree_1 = Phylo.read(path_1, tree_format)
    subtree_2 = Phylo.read(path_2, tree_format)

    list_1 = {i.name for i in subtree_1.get_terminals()}
    list_2 = {i.name for i in subtree_2.get_terminals()}

    size_list_1 = len(list_1)
    size_list_2 = len(list_2)

    intersection = list_1.intersection(list_2)

    return len(intersection) / max(size_list_1, size_list_2)

# %%
@python_app
def calcula_similaridade(max_rows: int, max_columns: int, matriz_subtree: list[list[Path]], tree_format: str) -> dict:
    dict_maf_database = {}

    for i in range(max_rows):
        for j in range(max_columns):
            for k in range(i + 1, max_rows):
                for l in range(max_columns):
                    g_maf = grade_maf(matriz_subtree[i][j], matriz_subtree[k][l], tree_format)

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

    # Falta validar as serquências
    analisar_sequencias_fasta("config.json")
    
    with open("config.json", "r", encoding="utf-8") as arquivo:
        parametros = json.load(arquivo)

    listas = dividir_lista(list(INPUT_SEQUENCES.iterdir()), NUCLEOS)    # Divide a lista em várias de acordo com o número de núcloes usados
    futuros = [alinhar_sequencias(parametros, lista) for lista in listas]   # Submete as tarefas ao Parsl
    wait(futuros)   # Espera todas terminarem

    listas = dividir_lista(list(SEQUENCIAS_ALINHADAS.iterdir()), NUCLEOS)    # Divide a lista em várias de acordo com o número de núcloes usados
    futuros = [constructor_tree(lista, parametros['tree_format']) for lista in listas]   # Submete as tarefas ao Parsl
    wait(futuros)   # Espera todas terminarem

    matriz_subtree = subarvores_possiveis(ARVORES_FILOGENETICAS, parametros['tree_format'])

    max_columns = max(len(row) for row in matriz_subtree)
    max_rows = len(matriz_subtree)

    matriz_subtree = preencher_matriz(matriz_subtree, None, max_columns)

    futures = calcula_similaridade(max_rows, max_columns, matriz_subtree, parametros['tree_format'])
    dict_maf_database = futures.result()

    # Salva o dicionário de similaridades
    with open(SIMILARIDADES / f"similaridades_{dt.datetime.now().strftime('%Y%m%d%H%M%S')}.json", "w", encoding="utf-8") as arquivo:
        json.dump(dict_maf_database, arquivo, ensure_ascii=False, indent=4)


    with open(PROVENANCE / "temp.json", "r", encoding="utf-8") as arquivo:
        par = json.load(arquivo)

    par['resultado'] = {"Inicio": inicio,
                        "Fim": time.perf_counter(),
                        "num_procs": NUCLEOS}
        
    par['host'] = system_summary()
        
    with open(PROVENANCE / f"dados_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.json", "w", encoding="utf-8") as arquivo:
        json.dump(par, arquivo, indent=4, ensure_ascii=False)
        
    os.remove(PROVENANCE / "temp.json")

# %%
clean_NoPipe()
clean_tmp()
clean_Trees()
clean_subtrees()

# %%



